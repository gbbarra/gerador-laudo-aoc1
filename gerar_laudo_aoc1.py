# -*- coding: utf-8 -*-

import pysam
import argparse
from collections import Counter
import datetime
import sys
import os

# --- BASE DE CONHECIMENTO (sem alterações) ---
VARIANTS = {
    'rs10156191': {
        'gene': 'AOC1', 'transcript': 'NM_001091.4', 'hgvsc': 'c.47C>T', 'hgvsp': 'p.Thr16Met', 'build': 'hg38', 'location': 'chr7:150856517', 'chrom': 'chr7', 'pos': 150856517, 'ref': 'C', 'alt': 'T',
        'genotypes': {
            'C/C': {'zygosity': 'Homozigoto Referência', 'result': 'padrão', 'interpretation': "O genótipo C/C (Thr/Thr) é o de referência (wild-type). Está associado à atividade normal da enzima DAO e ao risco basal para condições relacionadas à histamina."},
            'C/T': {'zygosity': 'Heterozigoto', 'result': 'moderadamente reduzida', 'interpretation': "O genótipo C/T (Thr/Met) está associado a uma redução moderada da atividade da DAO. Este resultado sugere uma predisposição aumentada a sintomas de intolerância à histamina, maior risco de hipersensibilidade a anti-inflamatórios não esteroides (AINEs) e enxaquecas, especialmente em mulheres."},
            'T/T': {'zygosity': 'Homozigoto Variante', 'result': 'reduzida', 'interpretation': "O genótipo T/T (Met/Met) está associado a uma redução ainda mais acentuada da atividade da DAO. Este resultado confere o maior risco para os sintomas mencionados, com forte predisposição à intolerância à histamina, manifestada por sintomas gastrointestinais e dores de cabeça."}
        }
    },
    'rs1049742': {
        'gene': 'AOC1', 'transcript': 'NM_001091.4', 'hgvsc': 'c.995C>T', 'hgvsp': 'p.Ser332Phe', 'build': 'hg38', 'location': 'chr7:150857465', 'chrom': 'chr7', 'pos': 150857465, 'ref': 'C', 'alt': 'T',
        'genotypes': {
            'C/C': {'zygosity': 'Homozigoto Referência', 'result': 'padrão', 'interpretation': "O genótipo C/C (Ser/Ser) é o de referência e está associado à atividade enzimática normal da DAO."},
            'C/T': {'zygosity': 'Heterozigoto', 'result': 'minimamente reduzida ou normal', 'interpretation': "O genótipo C/T (Ser/Phe) tem um efeito mínimo ou negligenciável na atividade da DAO. Geralmente, não está associado a um fenótipo clínico claro, mas pode contribuir para a intolerância à histamina apenas em combinação com outros fatores de risco."},
            'T/T': {'zygosity': 'Homozigoto Variante', 'result': 'levemente reduzida', 'interpretation': "O genótipo T/T (Phe/Phe) tem um efeito mínimo na atividade da DAO. Sendo raro e de baixo impacto, seu significado clínico é incerto e não está claramente associado a sintomas de intolerância à histamina de forma isolada."}
        }
    },
    'rs1049793': {
        'gene': 'AOC1', 'transcript': 'NM_001091.4', 'hgvsc': 'c.1933C>G', 'hgvsp': 'p.His645Asp', 'build': 'hg38', 'location': 'chr7:150860577', 'chrom': 'chr7', 'pos': 150860577, 'ref': 'C', 'alt': 'G',
        'genotypes': {
            'C/C': {'zygosity': 'Homozigoto Referência', 'result': 'padrão', 'interpretation': "O genótipo C/C (His/His) é o de referência e está associado à atividade enzimática normal da DAO."},
            'C/G': {'zygosity': 'Heterozigoto', 'result': 'reduzida (↓ ~34%)', 'interpretation': "O genótipo C/G (His/Asp) causa uma perda significativa da atividade da DAO (aprox. 34%). Este resultado indica um risco moderadamente aumentado para intolerância à histamina, com possível predisposição a sintomas gastrointestinais e cutâneos relacionados à histamina."},
            'G/G': {'zygosity': 'Homozigoto Variante', 'result': 'severamente reduzida (↓ ~49%)', 'interpretation': "O genótipo G/G (Asp/Asp) causa uma perda severa da atividade da DAO (aprox. 49%). Este resultado indica uma forte deficiência de DAO e um alto risco de intolerância à histamina, com predisposição a sintomas como distúrbios gastrointestinais, dores de cabeça e rubor facial."}
        }
    }
}

def genotype_position(bam_file, ref_file, chrom, pos, min_depth=10, min_base_quality=20, min_mapping_quality=20):
    """Realiza a genotipagem filtrando por Base Quality (BQ) e Mapping Quality (MQ)."""
    try:
        with pysam.AlignmentFile(bam_file, "rb") as samfile, \
             pysam.FastaFile(ref_file) as reffile:
            base_counts = Counter()
            for pileupcolumn in samfile.pileup(chrom, pos - 1, pos, fastafile=reffile, min_base_quality=min_base_quality, min_mapping_quality=min_mapping_quality):
                if pileupcolumn.pos == pos - 1:
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            base = pileupread.alignment.query_sequence[pileupread.query_position]
                            base_counts[base] += 1
    except FileNotFoundError as e:
        print(f"Erro: Arquivo não encontrado - {e}"); sys.exit(1)
    
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

def generate_html_report(results, sample_id, bam_filename, min_bq, min_mq):
    """Gera o laudo HTML objetivo com layout de resultado em duas linhas."""
    first_variant_info = next(iter(VARIANTS.values()))
    genome_build = first_variant_info['build']
    gene_name = first_variant_info['gene']
    transcript_id = first_variant_info['transcript']

    html_style = """
    <style>
        body { font-family: Arial, sans-serif; line-height: 1.6; color: #333; max-width: 900px; margin: 40px auto; }
        h1, h2 { color: #2c3e50; padding-bottom: 10px; }
        h1 { font-size: 24px; text-align: center; border-bottom: 2px solid #3498db; }
        h2 { font-size: 20px; margin-top: 40px; border-bottom: 1px solid #ecf0f1;}
        .patient-info { background-color: #f9f9f9; border-left: 5px solid #3498db; padding: 15px; margin-bottom: 15px; }
        .gene-info { 
            background-color: #f9f9f9; 
            padding: 15px 20px; 
            margin-bottom: 30px; 
            font-size: 16px; 
            text-align: left;
            line-height: 1.8;
            border-left: 5px solid #3498db;
        }
        .results-list { margin-top: 20px; }
        .variant-item { 
            background-color: #f8f9fa; 
            border: 1px solid #dee2e6; 
            padding: 15px; 
            margin-bottom: 12px; 
            border-radius: 4px;
            font-size: 15px;
            line-height: 1.7;
        }
        .variant-details {
            font-size: 14px;
            color: #555;
        }
        .interpretation-block { margin-bottom: 15px; padding-left: 10px; }
        strong { color: #000; }
        /* A classe footer foi removida para que os estilos sejam herdados do body */
        .note { font-style: italic; }
    </style>
    """

    # --- Pré-processamento dos resultados ---
    genotype_results_for_interpretation = {}
    any_variant_found = False
    for rsid, data in results.items():
        variant_info = VARIANTS[rsid]
        genotype_call = data['genotype']
        normalized_genotype = "/".join(sorted(genotype_call.split('/'))) if '/' in genotype_call else genotype_call
        
        genotype_results_for_interpretation[rsid] = {
            'full_hgv_name': f"{variant_info['transcript']}:{variant_info['hgvsc']} ({variant_info['hgvsp']})",
            'call': genotype_call, 'key': None
        }
        key = ""
        if "Baixa Cobertura" not in normalized_genotype:
            ref, alt = variant_info['ref'], variant_info['alt']
            if normalized_genotype == f"{ref}/{ref}": key = f"{ref}/{ref}"
            elif normalized_genotype == f"{ref}/{alt}": key = f"{ref}/{alt}"
            elif normalized_genotype == f"{alt}/{alt}": key = f"{alt}/{alt}"
        
        if key:
            genotype_results_for_interpretation[rsid]['key'] = key
            if 'Variante' in variant_info['genotypes'][key]['zygosity']:
                any_variant_found = True

    # --- Construção do HTML ---
    html_content = f"""
    <!DOCTYPE html><html lang="pt-br"><head><meta charset="UTF-8"><title>Laudo Genético - {sample_id}</title>{html_style}</head><body>
    <h1>Deficiência de diamina oxidase (DAO)</h1>
    <div class="patient-info">
        <strong>PACIENTE:</strong> {sample_id}<br>
        <strong>DATA DO LAUDO:</strong> {datetime.date.today().strftime('%d/%m/%Y')}<br>
        <strong>EXAME:</strong> Análise de variantes no gene AOC1 (DAO)<br>
        <strong>ARQUIVO ANALISADO:</strong> {os.path.basename(bam_filename)}
    </div>
    <div class="gene-info">
        <strong>Gene:</strong> {gene_name}<br>
        <strong>Transcrito:</strong> {transcript_id}<br>
        <strong>Referência genômica:</strong> {genome_build.upper()}
    </div>

    <h2>RESULTADO</h2>
    """

    if not any_variant_found:
        html_content += "<p>Nenhuma das variantes de risco analisadas foi detectada. O perfil genético do paciente é compatível com uma <strong>atividade normal</strong> da enzima DAO.</p>"
    else:
        html_content += "<p>Foi(ram) detectada(s) variante(s) que alteram a atividade da enzima DAO. O perfil genético do paciente sugere uma <strong>atividade enzimática reduzida</strong>.</p>"

    # --- Cria a lista de resultados com o novo layout de duas linhas ---
    results_list_html = '<div class="results-list">'
    
    for rsid, data in results.items():
        variant_info = VARIANTS[rsid]
        genotype_call = data['genotype']
        
        variant_hgv = f"<strong>{variant_info['hgvsc']} ({variant_info['hgvsp']})</strong>"
        location = f"{variant_info['location']} ({genome_build.upper()})"
        coverage_str = f"{data['coverage']}x"
        counts_str = ", ".join([f"{base} {count}x" for base, count in data['counts'].items()]) or "N/A"
        
        genotype_zygosity = "Indeterminado"
        activity = "Inconclusivo"
        
        if "Baixa Cobertura" in genotype_call:
            genotype_zygosity = f"<strong>{genotype_call}</strong>"
        else:
            key = genotype_results_for_interpretation[rsid]['key']
            if key:
                info = variant_info['genotypes'][key]
                genotype_zygosity = f"<strong>{key}</strong> ({info['zygosity']})"
                activity = f"<strong>{info['result']}</strong>"
            else:
                genotype_zygosity = f"<strong>{genotype_call}</strong> (Não canônico)"

        line1_html = f"{variant_hgv}  |  Genótipo: {genotype_zygosity}  |  Atividade Enzimática: {activity}"
        line2_html = f"ID: <strong>{rsid}</strong>  |  Posição: <strong>{location}</strong>  |  Cobertura (BQ≥{min_bq}): <strong>{coverage_str}</strong>  |  Leituras: <strong>{counts_str}</strong>"

        results_list_html += f"""
        <div class="variant-item">
            {line1_html}
            <div class="variant-details">{line2_html}</div>
        </div>
        """
    results_list_html += '</div>'
    html_content += results_list_html
    
    # Seção de Interpretação
    html_content += "<h2>INTERPRETAÇÃO DETALHADA</h2>"
    for rsid, result_data in genotype_results_for_interpretation.items():
        key = result_data.get('key')
        call = result_data['call']
        html_content += f"<div class='interpretation-block'><strong>• {result_data['full_hgv_name']}:</strong> "
        if "Baixa Cobertura" in call:
            html_content += "A análise foi <strong>inconclusiva</strong> devido à baixa cobertura de sequenciamento nesta posição, o que impede uma genotipagem confiável.</div>"
        elif key:
            html_content += f"{VARIANTS[rsid]['genotypes'][key]['interpretation']}</div>"
        else:
            html_content += f"O resultado para esta variante foi inconclusivo ({call}). O genótipo observado não corresponde aos alelos de referência/alternativo esperados.</div>"

    # --- SEÇÃO FINAL MODIFICADA (sem a classe .footer) ---
    html_content += f"""
        <h2>TÉCNICA APLICADA</h2>
        <p>A análise foi realizada a partir de DNA genômico extraído de forma automatizada (Maxwell®, Promega). O preparo da biblioteca de DNA foi conduzido com a tecnologia de captura em alvo "xGen™ DNA Library Prep EZ" e o painel de enriquecimento "xGen™ Exome Hyb Panel v2" (Integrated DNA Technologies - IDT), que abrange as regiões codificantes e intrônicas adjacentes (±10 pares de base) dos genes nucleares. O sequenciamento de nova geração (NGS) foi executado em uma plataforma MGI Tech DNABSEQ G400.</p>
        <p>O processamento primário dos dados brutos, incluindo o alinhamento das leituras ao genoma humano de referência (GRCh38/hg38) e a chamada de variantes, foi realizado através do pipeline Dragen Enrichment (v. 4.2.4).</p>
        <p>Para a elaboração deste laudo específico, as variantes de interesse no gene AOC1 foram submetidas a uma análise de genotipagem focada, executada por um pipeline próprio. Este processo consulta diretamente o arquivo de alinhamento (BAM) para verificar a cobertura e a identidade das bases em cada posição genômica alvo. A interpretação dos genótipos é baseada na base de conhecimento interna do script, compilada a partir das referências científicas listadas abaixo.</p>
        <p><strong>Limitações e Parâmetros de Qualidade:</strong> A análise se restringe apenas às variantes descritas neste laudo. A genotipagem depende da cobertura efetiva na posição de interesse (limite mínimo: 10x).</p>
        <p class="note"><strong>Nota sobre Qualidade (BQ):</strong> A genotipagem realizada considera apenas leituras com Qualidade de Base (Base Quality - BQ) igual ou superior a {min_bq}. A BQ é uma medida de confiança na identificação correta de cada base nucleotídica (A, C, G, T) pela plataforma de sequenciamento.</p>
        
        <h2>REFERÊNCIAS BIBLIOGRÁFICAS</h2>
        <p>1. Maintz L, et al. Allergy. 2011 Jul;66(7):893-902.<br>2. Ayuso P, et al. Pharmacogenet Genomics. 2007 Sep;17(9):687-93.<br>3. Agúndez JAG, et al. PLoS One. 2012;7(11):e47571.</p>
    </body></html>
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
    final_report = generate_html_report(
        results, args.sample_id, args.bam_file, args.min_base_quality, args.min_mapping_quality
    )
    
    with open(args.output_file, 'w', encoding='utf-8') as f:
        f.write(final_report)
    
    print(f"\nLaudo HTML salvo com sucesso em: {args.output_file}")
    print("Abra este arquivo em um navegador de internet para visualizar o laudo formatado.")


if __name__ == "__main__":
    main()