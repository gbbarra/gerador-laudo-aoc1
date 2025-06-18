# -*- coding: utf-8 -*-

import pysam
import argparse
from collections import Counter
import datetime
import sys

# --- CONFIGURAÇÃO DAS VARIANTES ---
# Coordenadas, build do genoma e transcrito atualizados conforme solicitação.
VARIANTS = {
    'rs10156191': {
        'name': 'p.Thr16Met',
        'build': 'hg38',
        'transcript': 'NM_001091.4',
        'location': 'chr7:150856517',
        'chrom': 'chr7',
        'pos': 150856517,
        'ref': 'C',
        'alt': 'T',
        'genotypes': {
            'C/C': {'zygosity': 'Homozigoto Referência', 'result': 'Atividade padrão'},
            'C/T': {'zygosity': 'Heterozigoto', 'result': 'Atividade moderadamente reduzida'},
            'T/T': {'zygosity': 'Homozigoto Variante', 'result': 'Atividade reduzida'},
        }
    },
    'rs1049742': {
        'name': 'p.Ser332Phe',
        'build': 'hg38',
        'transcript': 'NM_001091.4',
        'location': 'chr7:150857465',
        'chrom': 'chr7',
        'pos': 150857465,
        'ref': 'C',
        'alt': 'T',
        'genotypes': {
            'C/C': {'zygosity': 'Homozigoto Referência', 'result': 'Atividade padrão'},
            'C/T': {'zygosity': 'Heterozigoto', 'result': 'Atividade minimamente reduzida ou normal'},
            'T/T': {'zygosity': 'Homozigoto Variante', 'result': 'Atividade levemente reduzida'},
        }
    },
    'rs1049793': {
        'name': 'p.His645Asp',
        'build': 'hg38',
        'transcript': 'NM_001091.4',
        'location': 'chr7:150860577',
        'chrom': 'chr7',
        'pos': 150860577,
        'ref': 'C',
        'alt': 'G',
        'genotypes': {
            'C/C': {'zygosity': 'Homozigoto Referência', 'result': 'Atividade padrão'},
            'C/G': {'zygosity': 'Heterozigoto', 'result': 'Atividade reduzida (~34% de perda)'},
            'G/G': {'zygosity': 'Homozigoto Variante', 'result': 'Atividade severamente reduzida (~49% de perda)'},
        }
    }
}

def genotype_position(bam_file, ref_file, chrom, pos, min_depth=10, min_base_quality=20):
    """
    Realiza a genotipagem para uma única posição a partir de um arquivo BAM.
    """
    try:
        samfile = pysam.AlignmentFile(bam_file, "rb")
        if not samfile.has_index():
            print(f"Erro: Arquivo de índice (.bai) não encontrado para {bam_file}. Execute 'samtools index {bam_file}'.")
            sys.exit(1)
        reffile = pysam.FastaFile(ref_file)
    except FileNotFoundError as e:
        print(f"Erro: Arquivo não encontrado - {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Erro ao abrir os arquivos: {e}. Verifique se os nomes dos cromossomos são consistentes.")
        sys.exit(1)

    base_counts = Counter()
    coverage = 0
    try:
        for pileupcolumn in samfile.pileup(chrom, pos - 1, pos, fastafile=reffile, min_base_quality=min_base_quality):
            if pileupcolumn.pos == pos - 1:
                coverage = pileupcolumn.nsegments
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        base_counts[base] += 1
    except ValueError as e:
        print(f"\nErro Crítico: O cromossomo '{chrom}' não foi encontrado no arquivo de referência '{ref_file}'.")
        print("Verifique se o seu arquivo BAM e o arquivo de referência usam o mesmo formato de cromossomos (ex: 'chr7' vs '7').")
        samfile.close()
        reffile.close()
        sys.exit(1)

    samfile.close()
    reffile.close()

    if coverage < min_depth:
        return "Baixa Cobertura", coverage, base_counts
    if not base_counts:
        return "Sem Leitura", coverage, base_counts

    total_reads = sum(base_counts.values())
    top_alleles = [item[0] for item in base_counts.most_common(2)]
    
    if len(top_alleles) == 1:
        return f"{top_alleles[0]}/{top_alleles[0]}", coverage, base_counts
    else:
        allele1, allele2 = top_alleles[0], top_alleles[1]
        freq_allele2 = base_counts[allele2] / total_reads
        if 0.20 < freq_allele2 < 0.80:
             return "/".join(sorted(top_alleles)), coverage, base_counts
        else:
            return f"{allele1}/{allele1}", coverage, base_counts


def generate_report(results, sample_id="N/A"):
    """
    Gera o laudo final com base nos resultados da genotipagem.
    """
    # Pega o build e o transcrito da primeira variante para usar no texto do laudo
    first_variant_info = next(iter(VARIANTS.values()))
    genome_build = first_variant_info['build']
    transcript_id = first_variant_info['transcript']

    report = f"""
================================================================================
                                LAUDO GENÉTICO
================================================================================

PACIENTE: {sample_id}
DATA DO LAUDO: {datetime.date.today().strftime('%d/%m/%Y')}
EXAME: Análise de variantes no gene AOC1 (DAO) por Sequenciamento de Nova Geração

--------------------------------------------------------------------------------
INTRODUÇÃO: GENE AOC1 (DAO) E METABOLISMO DA HISTAMINA
--------------------------------------------------------------------------------
O gene AOC1 (diamina oxidase) codifica a principal enzima responsável pela degradação da histamina no trato digestivo. Variantes genéticas neste gene podem reduzir a atividade da enzima, levando a um acúmulo de histamina e a um quadro clínico conhecido como Intolerância à Histamina (HIT), que pode incluir sintomas como enxaqueca, distúrbios gastrointestinais e reações cutâneas.

--------------------------------------------------------------------------------
RESULTADOS
--------------------------------------------------------------------------------
"""
    table_header = f"{'Variante Analisada':<20} {'Id. SNV':<15} {'Posição Genômica ({genome_build.upper()})':<22} {'Genótipo':<10} {'Zigosidade':<25} {'Resultado (Atividade Enzimática)'}"
    report += table_header + "\n"
    report += "-" * len(table_header) + "\n"

    found_variants_info = []
    for rsid, data in results.items():
        variant_info = VARIANTS[rsid]
        location_str = variant_info['location']
        genotype_call = data['genotype']
        
        alleles = sorted(genotype_call.split('/')) if '/' in genotype_call else [genotype_call]
        normalized_genotype = "/".join(alleles)

        if "Baixa Cobertura" in normalized_genotype or "Sem Leitura" in normalized_genotype:
            zygosity = "Indeterminado"
            result_text = genotype_call
            genotype_display = "N/A"
        else:
            genotype_key = f"{variant_info['ref']}/{variant_info['alt']}"
            if normalized_genotype == f"{variant_info['ref']}/{variant_info['ref']}":
                genotype_key = f"{variant_info['ref']}/{variant_info['ref']}"
            elif normalized_genotype == f"{variant_info['alt']}/{variant_info['alt']}":
                 genotype_key = f"{variant_info['alt']}/{variant_info['alt']}"
            
            genotype_info = variant_info['genotypes'].get(genotype_key, {'zygosity': 'Desconhecido', 'result': 'Genótipo não canônico'})
            
            zygosity = genotype_info['zygosity']
            result_text = genotype_info['result']
            genotype_display = genotype_key

        report += f"{variant_info['name']:<20} {rsid:<15} {location_str:<22} {genotype_display:<10} {zygosity:<25} {result_text}\n"

        if "Variante" in zygosity:
            found_variants_info.append({'rsid': rsid, 'name': variant_info['name'], 'zygosity': zygosity.lower().replace(' variante', '')})

    report += """
--------------------------------------------------------------------------------
INTERPRETAÇÃO DO RESULTADO
--------------------------------------------------------------------------------
"""
    if not found_variants_info:
        report += f"NÃO FOI DETECTADA a presença de nenhuma das variantes analisadas (p.Thr16Met, p.Ser332Phe, p.His645Asp) no gene AOC1 (transcrito {transcript_id}). O genótipo do paciente é o de referência ('wild-type'), o que é compatível com uma atividade enzimática padrão da DAO.\n"
    else:
        variant_summary_parts = [f"a variante {v['name']} ({v['rsid']}) em {v['zygosity']}" for v in found_variants_info]
        report += f"Foi detectada a presença de: {', '.join(variant_summary_parts)}.\n\n"
        if any(v['rsid'] == 'rs10156191' for v in found_variants_info):
            report += "- **p.Thr16Met (rs10156191)**: Associada a uma redução moderada da atividade da DAO. Portadores do alelo variante podem ter maior risco de intolerância à histamina, enxaquecas e hipersensibilidade a AINEs.\n"
        if any(v['rsid'] == 'rs1049742' for v in found_variants_info):
            report += "- **p.Ser332Phe (rs1049742)**: O impacto desta variante na atividade da DAO é considerado mínimo. Geralmente, não está associada a um fenótipo clínico claro.\n"
        if any(v['rsid'] == 'rs1049793' for v in found_variants_info):
             report += "- **p.His645Asp (rs1049793)**: Associada a uma redução substancial na atividade da DAO. Portadores do alelo variante têm um risco elevado de deficiência de DAO e sintomas associados.\n"

    report += """
Este laudo deve ser interpretado por um especialista dentro do contexto clínico e familiar do paciente.

--------------------------------------------------------------------------------
TÉCNICA APLICADA E LIMITAÇÕES
--------------------------------------------------------------------------------
O DNA genômico foi analisado por Sequenciamento de Nova Geração (NGS). As sequências foram alinhadas ao genoma humano de referência ({genome_build.upper()}) e as variantes de interesse no gene AOC1, baseadas no transcrito de referência {transcript_id}, foram genotipadas. A análise se restringe às variantes descritas. A genotipagem depende da cobertura na posição de interesse (limite: 10x). Variantes estruturais complexas não são detectadas.

--------------------------------------------------------------------------------
REFERÊNCIAS BIBLIOGRÁFICAS
--------------------------------------------------------------------------------
1. Maintz L, et al. Allergy. 2011 Jul;66(7):893-902.
2. Ayuso P, et al. Pharmacogenet Genomics. 2007 Sep;17(9):687-93.
3. Agúndez JAG, et al. PLoS One. 2012;7(11):e47571.

================================================================================
"""
    return report

def main():
    parser = argparse.ArgumentParser(description="Gera um laudo para variantes do gene AOC1 a partir de um arquivo BAM.")
    parser.add_argument("bam_file", help="Caminho para o arquivo BAM do paciente.")
    parser.add_argument("ref_file", help="Caminho para o arquivo FASTA do genoma de referência.")
    parser.add_argument("--sample_id", help="ID do paciente/amostra para o laudo.", default="Amostra Anônima")
    parser.add_argument("--output_file", help="Nome do arquivo de texto para salvar o laudo.", default=None)
    
    args = parser.parse_args()

    results = {}
    print("Analisando variantes no gene AOC1...")
    for rsid, info in VARIANTS.items():
        print(f"  - Genotipando {info['name']} ({rsid}) em {info['location']} (Transcrito: {info['transcript']})...")
        genotype, coverage, counts = genotype_position(args.bam_file, args.ref_file, info['chrom'], info['pos'])
        results[rsid] = {'genotype': genotype, 'coverage': coverage, 'counts': dict(counts)}
        print(f"    -> Resultado: Genótipo={genotype}, Cobertura={coverage}x, Contagens de Bases={dict(counts)}")

    print("\nGerando laudo...")
    final_report = generate_report(results, args.sample_id)
    
    if args.output_file:
        with open(args.output_file, 'w', encoding='utf-8') as f:
            f.write(final_report)
        print(f"\nLaudo salvo com sucesso em: {args.output_file}")
    else:
        print("\n--- INÍCIO DO LAUDO ---")
        print(final_report)
        print("--- FIM DO LAUDO ---")

if __name__ == "__main__":
    main()