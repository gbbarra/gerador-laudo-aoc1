# -*- coding: utf-8 -*-

import pysam
import argparse
from collections import Counter
import datetime
import sys

# --- BASE DE CONHECIMENTO (sem alterações) ---
VARIANTS = {
    'rs10156191': {
        'gene': 'AOC1', 'transcript': 'NM_001091.4', 'hgvsc': 'c.47C>T', 'hgvsp': 'p.Thr16Met', 'build': 'hg38', 'location': 'chr7:150856517', 'chrom': 'chr7', 'pos': 150856517, 'ref': 'C', 'alt': 'T',
        'genotypes': {
            'C/C': {'zygosity': 'Homozigoto Referência', 'result': 'Atividade padrão', 'interpretation': "O genótipo C/C (Thr/Thr) é o de referência (wild-type). Está associado à atividade normal da enzima DAO e ao risco basal para condições relacionadas à histamina."},
            'C/T': {'zygosity': 'Heterozigoto', 'result': 'Atividade moderadamente reduzida', 'interpretation': "O genótipo C/T (Thr/Met) está associado a uma redução moderada da atividade da DAO. Este resultado sugere uma predisposição aumentada a sintomas de intolerância à histamina, maior risco de hipersensibilidade a anti-inflamatórios não esteroides (AINEs) e enxaquecas, especialmente em mulheres."},
            'T/T': {'zygosity': 'Homozigoto Variante', 'result': 'Atividade reduzida', 'interpretation': "O genótipo T/T (Met/Met) está associado a uma redução ainda mais acentuada da atividade da DAO. Este resultado confere o maior risco para os sintomas mencionados, com forte predisposição à intolerância à histamina, manifestada por sintomas gastrointestinais e dores de cabeça."}
        }
    },
    'rs1049742': {
        'gene': 'AOC1', 'transcript': 'NM_001091.4', 'hgvsc': 'c.995C>T', 'hgvsp': 'p.Ser332Phe', 'build': 'hg38', 'location': 'chr7:150857465', 'chrom': 'chr7', 'pos': 150857465, 'ref': 'C', 'alt': 'T',
        'genotypes': {
            'C/C': {'zygosity': 'Homozigoto Referência', 'result': 'Atividade padrão', 'interpretation': "O genótipo C/C (Ser/Ser) é o de referência e está associado à atividade enzimática normal da DAO."},
            'C/T': {'zygosity': 'Heterozigoto', 'result': 'Atividade minimamente reduzida ou normal', 'interpretation': "O genótipo C/T (Ser/Phe) tem um efeito mínimo ou negligenciável na atividade da DAO. Geralmente, não está associado a um fenótipo clínico claro, mas pode contribuir para a intolerância à histamina apenas em combinação com outros fatores de risco."},
            'T/T': {'zygosity': 'Homozigoto Variante', 'result': 'Atividade levemente reduzida', 'interpretation': "O genótipo T/T (Phe/Phe) tem um efeito mínimo na atividade da DAO. Sendo raro e de baixo impacto, seu significado clínico é incerto e não está claramente associado a sintomas de intolerância à histamina de forma isolada."}
        }
    },
    'rs1049793': {
        'gene': 'AOC1', 'transcript': 'NM_001091.4', 'hgvsc': 'c.1933C>G', 'hgvsp': 'p.His645Asp', 'build': 'hg38', 'location': 'chr7:150860577', 'chrom': 'chr7', 'pos': 150860577, 'ref': 'C', 'alt': 'G',
        'genotypes': {
            'C/C': {'zygosity': 'Homozigoto Referência', 'result': 'Atividade padrão', 'interpretation': "O genótipo C/C (His/His) é o de referência e está associado à atividade enzimática normal da DAO."},
            'C/G': {'zygosity': 'Heterozigoto', 'result': 'Atividade reduzida (~34% de perda)', 'interpretation': "O genótipo C/G (His/Asp) causa uma perda significativa da atividade da DAO (aprox. 34%). Este resultado indica um risco moderadamente aumentado para intolerância à histamina, com possível predisposição a sintomas gastrointestinais e cutâneos relacionados à histamina."},
            'G/G': {'zygosity': 'Homozigoto Variante', 'result': 'Atividade severamente reduzida (~49% de perda)', 'interpretation': "O genótipo G/G (Asp/Asp) causa uma perda severa da atividade da DAO (aprox. 49%). Este resultado indica uma forte deficiência de DAO e um alto risco de intolerância à histamina, com predisposição a sintomas como distúrbios gastrointestinais, dores de cabeça e rubor facial."}
        }
    }
}

def genotype_position(bam_file, ref_file, chrom, pos, min_depth=10, min_base_quality=20, min_mapping_quality=20):
    """Realiza a genotipagem filtrando por Base Quality (BQ) e Mapping Quality (MQ)."""
    # (O código desta função permanece o mesmo)
    try:
        samfile = pysam.AlignmentFile(bam_file, "rb")
        reffile = pysam.FastaFile(ref_file)
    except FileNotFoundError as e:
        print(f"Erro: Arquivo não encontrado - {e}"); sys.exit(1)
    base_counts = Counter()
    for pileupcolumn in samfile.pileup(chrom, pos - 1, pos, fastafile=reffile, min_base_quality=min_base_quality, min_mapping_quality=min_mapping_quality):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base_counts[pileupread.alignment.query_sequence[pileupread.query_position]] += 1
    samfile.close(); reffile.close()
    effective_coverage = sum(base_counts.values())
    if effective_coverage < min_depth:
        return "Baixa Cobertura", effective_coverage, base_counts
    top_alleles = [item[0] for item in base_counts.most_common(2)]
    if len(top_alleles) == 1:
        return f"{top_alleles[0]}/{top_alleles[0]}", effective_coverage, base_counts
    else:
        allele1, allele2 = top_alleles[0], top_alleles[1]
        freq_allele2 = base_counts[allele2] / effective_coverage
        if 0.20 < freq_allele2 < 0.80:
            return "/".join(sorted(top_alleles)), effective_coverage, base_counts
        else:
            return f"{allele1}/{allele1}", effective_coverage, base_counts

def generate_html_report(results, sample_id, min_bq, min_mq):
    """Gera o laudo HTML final, com a explicação sobre BQ."""
    first_variant_info = next(iter(VARIANTS.values()))
    genome_build = first_variant_info['build']
    transcript_id = first_variant_info['transcript']

    html_style = """
    <style>
        body { font-family: Arial, sans-serif; line-height: 1.6; color: #333; max-width: 800px; margin: 40px auto; }
        h1, h2, h3 { color: #2c3e50; padding-bottom: 10px; }
        h1 { font-size: 24px; text-align: center; border-bottom: 2px solid #3498db; }
        h2 { font-size: 20px; margin-top: 40px; border-bottom: 1px solid #ecf0f1;}
        h3 { font-size: 18px; margin-top: 30px; border-bottom: none; color: #34495e; }
        .patient-info { background-color: #f9f9f9; border-left: 5px solid #3498db; padding: 15px; margin-bottom: 30px; }
        .result-block { margin-bottom: 30px; }
        .result-item { padding-left: 10px; }
        .interpretation-block { margin-bottom: 15px; }
        strong { color: #000; }
        .footer { font-size: 12px; color: #777; margin-top: 50px; text-align: left; border-top: 1px solid #ecf0f1; padding-top: 20px;}
        .note { font-style: italic; }
    </style>
    """

    result_blocks_html = ""
    genotype_results_for_interpretation = {}
    for rsid, data in results.items():
        variant_info = VARIANTS[rsid]
        full_hgv_name = f"<strong>{variant_info['gene']}</strong>:{variant_info['transcript']}:{variant_info['hgvsc']} ({variant_info['hgvsp']})"
        genotype_results_for_interpretation[rsid] = {'full_hgv_name': full_hgv_name.replace('<strong>', '').replace('</strong>', '')}

        genotype_call = data['genotype']
        normalized_genotype = "/".join(sorted(genotype_call.split('/'))) if '/' in genotype_call else genotype_call
        genotype_results_for_interpretation[rsid]['call'] = normalized_genotype
        
        block = f"<div class='result-block'><h3>{full_hgv_name}</h3>"
        block += "<div class='result-item'>"
        block += f"ID (dbSNP): {rsid}<br>"
        block += f"Posição Genômica ({genome_build.upper()}): {variant_info['location']}<br>"
        
        if "Baixa Cobertura" in normalized_genotype:
            block += f"Genótipo Identificado: <strong>INDETERMINADO ({genotype_call})</strong><br>"
            block += f"Cobertura Efetiva: <strong>{data['coverage']}x</strong> (leituras com BQ ≥ {min_bq})"
        else:
            ref, alt = variant_info['ref'], variant_info['alt']
            key = ""
            if normalized_genotype == f"{ref}/{ref}": key = f"{ref}/{ref}"
            elif normalized_genotype == f"{ref}/{alt}": key = f"{ref}/{alt}"
            elif normalized_genotype == f"{alt}/{alt}": key = f"{alt}/{alt}"

            if key:
                info = variant_info['genotypes'][key]
                genotype_results_for_interpretation[rsid]['key'] = key
                counts_str = ", ".join([f"{base} {count}x" for base, count in data['counts'].items()])

                block += f"Genótipo Identificado: <strong>{key} ({info['zygosity']})</strong><br>"
                block += f"Cobertura Efetiva: <strong>{data['coverage']}x</strong> (leituras com BQ ≥ {min_bq}) | Leituras: <strong>{counts_str}</strong><br>"
                block += f"Resultado: <strong>{info['result']}</strong>"
            else:
                block += f"Genótipo Identificado: {normalized_genotype} (Genótipo não canônico)"
        block += "</div></div>"
        result_blocks_html += block
    
    any_variant_found = any('Variante' in VARIANTS[rsid]['genotypes'].get(res.get('key', ''), {}).get('zygosity', '') for rsid, res in genotype_results_for_interpretation.items())
    
    html_content = f"""
    <!DOCTYPE html><html lang="pt-br"><head><meta charset="UTF-8"><title>Laudo Genético - {sample_id}</title>{html_style}</head><body>
    <h1>LAUDO GENÉTICO</h1>
    <div class="patient-info"><strong>PACIENTE:</strong> {sample_id}<br><strong>DATA DO LAUDO:</strong> {datetime.date.today().strftime('%d/%m/%Y')}<br><strong>EXAME:</strong> Análise de variantes no gene AOC1 (DAO)</div>
    <h2>INTRODUÇÃO</h2><p>O gene AOC1 (diamina oxidase) codifica a principal enzima responsável pela degradação da histamina no trato digestivo. Variantes genéticas neste gene podem reduzir a atividade da enzima, levando a um acúmulo de histamina e a um quadro clínico conhecido como Intolerância à Histamina (HIT), que pode incluir sintomas como enxaqueca, distúrbios gastrointestinais e reações cutâneas.</p>
    <h2>RESULTADOS DA ANÁLISE</h2>{result_blocks_html}
    <h2>INTERPRETAÇÃO DETALHADA</h2>
    """
    if not any_variant_found: html_content += "<p><strong>RESUMO GERAL:</strong> Nenhuma das variantes de risco analisadas foi detectada. O perfil genético do paciente é compatível com uma <strong>atividade normal</strong> da enzima DAO.</p>"
    else: html_content += "<p><strong>RESUMO GERAL:</strong> Foi(ram) detectada(s) variante(s) que alteram a atividade da enzima DAO. O perfil genético do paciente sugere uma <strong>atividade enzimática reduzida</strong>. A análise detalhada abaixo descreve o impacto específico de cada genótipo.</p>"
    
    html_content += "<h3>ANÁLISE POR VARIANTE</h3>"
    for rsid, result_data in genotype_results_for_interpretation.items():
        key = result_data.get('key')
        html_content += f"<div class='interpretation-block'><strong>• {result_data['full_hgv_name']}:</strong> "
        if key: html_content += f"{VARIANTS[rsid]['genotypes'][key]['interpretation']}</div>"
        else: html_content += f"O resultado para esta variante foi inconclusivo ({result_data['call']}).</div>"
        
    html_content += f"""
    <div class="footer">
        <h2>TÉCNICA APLICADA E LIMITAÇÕES</h2>
        <p>O DNA genômico foi analisado por Sequenciamento de Nova Geração (NGS). As sequências foram alinhadas ao genoma humano de referência ({genome_build.upper()}) e as variantes de interesse no gene AOC1, baseadas no transcrito de referência {transcript_id}, foram genotipadas. A análise se restringe às variantes descritas. A genotipagem depende da cobertura efetiva na posição de interesse (limite: 10x).</p>
        <p class="note"><strong>Nota sobre Qualidade (BQ):</strong> A genotipagem considera apenas leituras com Qualidade de Base (Base Quality - BQ) igual ou superior a {min_bq}. A BQ é uma medida de confiança na identificação correta de cada base nucleotídica (A, C, G, T) pela plataforma de sequenciamento.</p>
        
        <h2>REFERÊNCIAS BIBLIOGRÁFICAS</h2>
        <p>1. Maintz L, et al. Allergy. 2011 Jul;66(7):893-902.<br>2. Ayuso P, et al. Pharmacogenet Genomics. 2007 Sep;17(9):687-93.<br>3. Agúndez JAG, et al. PLoS One. 2012;7(11):e47571.</p>
    </div></body></html>
    """
    return html_content

def main():
    parser = argparse.ArgumentParser(description="Gera um laudo para variantes do gene AOC1 a partir de um arquivo BAM.")
    parser.add_argument("bam_file", help="Caminho para o arquivo BAM do paciente.")
    parser.add_argument("ref_file", help="Caminho para o arquivo FASTA do genoma de referência.")
    parser.add_argument("--sample_id", help="ID do paciente/amostra para o laudo.", default="Amostra Anônima")
    parser.add_argument("--output_file", help="Nome do arquivo HTML para salvar o laudo (ex: laudo.html).", default="laudo_aoc1.html")
    parser.add_argument("--min_base_quality", type=int, default=20, help="Qualidade mínima de base (BQ) para considerar uma leitura.")
    parser.add_argument("--min_mapping_quality", type=int, default=20, help="Qualidade mínima de mapeamento (MQ) para considerar uma leitura.")
    
    args = parser.parse_args()
    if not args.output_file.lower().endswith('.html'): args.output_file += '.html'
    
    results = {}
    print("Analisando variantes no gene AOC1...")
    for rsid, info in VARIANTS.items():
        print(f"  - Genotipando {info['gene']}:{info['hgvsp']} ({rsid}) com BQ>={args.min_base_quality}, MQ>={args.min_mapping_quality}...")
        genotype, effective_coverage, counts = genotype_position(
            args.bam_file, args.ref_file, info['chrom'], info['pos'], 
            min_base_quality=args.min_base_quality, 
            min_mapping_quality=args.min_mapping_quality
        )
        results[rsid] = {'genotype': genotype, 'coverage': effective_coverage, 'counts': dict(counts)}
        print(f"    -> Resultado: Genótipo={genotype}, Cobertura Efetiva={effective_coverage}x, Contagens={dict(counts)}")
    
    print("\nGerando laudo em HTML...")
    final_report = generate_html_report(results, args.sample_id, args.min_base_quality, args.min_mapping_quality)
    with open(args.output_file, 'w', encoding='utf-8') as f: f.write(final_report)
    
    print(f"\nLaudo HTML salvo com sucesso em: {args.output_file}")
    print("Abra este arquivo em um navegador de internet para visualizar o laudo formatado.")


if __name__ == "__main__":
    main()