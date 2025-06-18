# -*- coding: utf-8 -*-

import pysam
import argparse
from collections import Counter
import datetime
import sys

# --- CONFIGURAÇÃO DAS VARIANTES COM INTERPRETAÇÕES DETALHADAS POR GENÓTIPO ---
# A base de conhecimento permanece a mesma.
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
            'C/C': {
                'zygosity': 'Homozigoto Referência',
                'result': 'Atividade padrão',
                'interpretation': "O genótipo C/C (Thr/Thr) é o de referência (wild-type). Está associado à atividade normal da enzima DAO e ao risco basal para condições relacionadas à histamina."
            },
            'C/T': {
                'zygosity': 'Heterozigoto',
                'result': 'Atividade moderadamente reduzida',
                'interpretation': "O genótipo C/T (Thr/Met) está associado a uma redução moderada da atividade da DAO. Este resultado sugere uma predisposição aumentada a sintomas de intolerância à histamina, maior risco de hipersensibilidade a anti-inflamatórios não esteroides (AINEs) e enxaquecas, especialmente em mulheres."
            },
            'T/T': {
                'zygosity': 'Homozigoto Variante',
                'result': 'Atividade reduzida',
                'interpretation': "O genótipo T/T (Met/Met) está associado a uma redução ainda mais acentuada da atividade da DAO. Este resultado confere o maior risco para os sintomas mencionados, com forte predisposição à intolerância à histamina, manifestada por sintomas gastrointestinais e dores de cabeça."
            }
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
            'C/C': {
                'zygosity': 'Homozigoto Referência',
                'result': 'Atividade padrão',
                'interpretation': "O genótipo C/C (Ser/Ser) é o de referência e está associado à atividade enzimática normal da DAO."
            },
            'C/T': {
                'zygosity': 'Heterozigoto',
                'result': 'Atividade minimamente reduzida ou normal',
                'interpretation': "O genótipo C/T (Ser/Phe) tem um efeito mínimo ou negligenciável na atividade da DAO. Geralmente, não está associado a um fenótipo clínico claro, mas pode contribuir para a intolerância à histamina apenas em combinação com outros fatores de risco."
            },
            'T/T': {
                'zygosity': 'Homozigoto Variante',
                'result': 'Atividade levemente reduzida',
                'interpretation': "O genótipo T/T (Phe/Phe) tem um efeito mínimo na atividade da DAO. Sendo raro e de baixo impacto, seu significado clínico é incerto e não está claramente associado a sintomas de intolerância à histamina de forma isolada."
            }
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
            'C/C': {
                'zygosity': 'Homozigoto Referência',
                'result': 'Atividade padrão',
                'interpretation': "O genótipo C/C (His/His) é o de referência e está associado à atividade enzimática normal da DAO."
            },
            'C/G': {
                'zygosity': 'Heterozigoto',
                'result': 'Atividade reduzida (~34% de perda)',
                'interpretation': "O genótipo C/G (His/Asp) causa uma perda significativa da atividade da DAO (aprox. 34%). Este resultado indica um risco moderadamente aumentado para intolerância à histamina, com possível predisposição a sintomas gastrointestinais e cutâneos relacionados à histamina."
            },
            'G/G': {
                'zygosity': 'Homozigoto Variante',
                'result': 'Atividade severamente reduzida (~49% de perda)',
                'interpretation': "O genótipo G/G (Asp/Asp) causa uma perda severa da atividade da DAO (aprox. 49%). Este resultado indica uma forte deficiência de DAO e um alto risco de intolerância à histamina, com predisposição a sintomas como distúrbios gastrointestinais, dores de cabeça e rubor facial."
            }
        }
    }
}

def genotype_position(bam_file, ref_file, chrom, pos, min_depth=10, min_base_quality=20):
    """Realiza a genotipagem para uma única posição a partir de um arquivo BAM."""
    # (O código desta função permanece o mesmo)
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
    """Gera o laudo com o novo estilo visual: sem linhas, usando negrito e maiúsculas."""
    first_variant_info = next(iter(VARIANTS.values()))
    genome_build = first_variant_info['build']
    transcript_id = first_variant_info['transcript']

    # --- Seção de Resultados em formato de lista ---
    result_blocks = []
    genotype_results_for_interpretation = {}
    for rsid, data in results.items():
        variant_info = VARIANTS[rsid]
        location_str = variant_info['location']
        genotype_call = data['genotype']
        
        alleles = sorted(genotype_call.split('/')) if '/' in genotype_call else [genotype_call]
        normalized_genotype = "/".join(alleles)
        genotype_results_for_interpretation[rsid] = {'call': normalized_genotype}

        block = f"**{variant_info['name'].upper()} ({rsid})**\n"
        block += f"  - Posição Genômica ({genome_build.upper()}): {location_str}\n"

        if "Baixa Cobertura" in normalized_genotype or "Sem Leitura" in normalized_genotype:
            block += f"  - Genótipo Identificado: INDETERMINADO ({genotype_call})"
        else:
            ref_allele, alt_allele = variant_info['ref'], variant_info['alt']
            genotype_key = ""
            if normalized_genotype == f"{ref_allele}/{ref_allele}": genotype_key = f"{ref_allele}/{ref_allele}"
            elif normalized_genotype == f"{ref_allele}/{alt_allele}": genotype_key = f"{ref_allele}/{alt_allele}"
            elif normalized_genotype == f"{alt_allele}/{alt_allele}": genotype_key = f"{alt_allele}/{alt_allele}"

            if genotype_key:
                genotype_info = variant_info['genotypes'][genotype_key]
                genotype_results_for_interpretation[rsid]['key'] = genotype_key
                block += f"  - Genótipo Identificado: **{genotype_key} ({genotype_info['zygosity']})**\n"
                block += f"  - Resultado: **{genotype_info['result']}**"
            else:
                block += f"  - Genótipo Identificado: {normalized_genotype} (Genótipo não canônico)"
        
        result_blocks.append(block)

    # --- Lógica do Resumo Geral ---
    any_variant_found = any('Variante' in VARIANTS[rsid]['genotypes'].get(res.get('key', ''), {}).get('zygosity', '') for rsid, res in genotype_results_for_interpretation.items())

    # --- Construção do Laudo Final ---
    report = f"""
**LAUDO GENÉTICO**

**INFORMAÇÕES DO PACIENTE**
PACIENTE: {sample_id}
DATA DO LAUDO: {datetime.date.today().strftime('%d/%m/%Y')}
EXAME: Análise de variantes no gene AOC1 (DAO)

**INTRODUÇÃO**
O gene AOC1 (diamina oxidase) codifica a principal enzima responsável pela degradação da histamina no trato digestivo. Variantes genéticas neste gene podem reduzir a atividade da enzima, levando a um acúmulo de histamina e a um quadro clínico conhecido como Intolerância à Histamina (HIT), que pode incluir sintomas como enxaqueca, distúrbios gastrointestinais e reações cutâneas.

**RESULTADOS DA ANÁLISE**
{'\n\n'.join(result_blocks)}

**INTERPRETAÇÃO DETALHADA**
"""
    if not any_variant_found:
        report += f"**RESUMO GERAL:** Nenhuma das variantes de risco analisadas foi detectada. O perfil genético do paciente é compatível com uma **atividade normal** da enzima DAO.\n\n"
    else:
        report += f"**RESUMO GERAL:** Foi(ram) detectada(s) variante(s) que alteram a atividade da enzima DAO. O perfil genético do paciente sugere uma **atividade enzimática reduzida**. A análise detalhada abaixo descreve o impacto específico de cada genótipo.\n\n"

    report += "**ANÁLISE POR VARIANTE:**\n"
    for rsid, result_data in genotype_results_for_interpretation.items():
        variant_info = VARIANTS[rsid]
        report += f"• **{variant_info['name']} ({rsid}):** "
        genotype_key = result_data.get('key')
        if genotype_key:
            report += f"{variant_info['genotypes'][genotype_key]['interpretation']}\n"
        else:
            report += f"O resultado para esta variante foi inconclusivo ({result_data['call']}).\n"

    report += f"""
**TÉCNICA APLICADA E LIMITAÇÕES**
O DNA genômico foi analisado por Sequenciamento de Nova Geração (NGS). As sequências foram alinhadas ao genoma humano de referência ({genome_build.upper()}) e as variantes de interesse no gene AOC1, baseadas no transcrito de referência {transcript_id}, foram genotipadas. A análise se restringe às variantes descritas. A genotipagem depende da cobertura na posição de interesse (limite: 10x). Variantes estruturais complexas não são detectadas.

**REFERÊNCIAS BIBLIOGRÁFICAS**
1. Maintz L, et al. Allergy. 2011 Jul;66(7):893-902.
2. Ayuso P, et al. Pharmacogenet Genomics. 2007 Sep;17(9):687-93.
3. Agúndez JAG, et al. PLoS One. 2012;7(11):e47571.
"""
    return report

def main():
    # (O código desta função permanece o mesmo)
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