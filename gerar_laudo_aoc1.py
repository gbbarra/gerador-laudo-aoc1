# -*- coding: utf-8 -*-

import pysam
import argparse
from collections import Counter
import datetime
import sys
import os

# --- BASE DE CONHECIMENTO (Ajustada para separar resultado funcional da interpretação clínica) ---
VARIANTS = {
    'rs10156191': {
        'gene': 'AOC1', 'transcript': 'NM_001091.4', 'hgvsc': 'c.47C>T', 'hgvsp': 'p.Thr16Met', 'build': 'hg38', 'location': 'chr7:150856517', 'chrom': 'chr7', 'pos': 150856517, 'ref': 'C', 'alt': 'T',
        'genotypes': {
            'C/C': {'zygosity': 'Homozigoto Referência', 'result': 'padrão', 
                    'interpretation': "O genótipo C/C (Thr/Thr) é o de referência (wild-type) e está associado à atividade normal da enzima DAO."},
            'C/T': {'zygosity': 'Heterozigoto', 'result': 'moderadamente reduzida', 
                    'interpretation': "O genótipo C/T (Thr/Met) está associado a uma redução moderada da atividade da DAO."},
            'T/T': {'zygosity': 'Homozigoto Variante', 'result': 'reduzida', 
                    'interpretation': "O genótipo T/T (Met/Met) está associado a uma redução na atividade da DAO."}
        }
    },
    'rs1049742': {
        'gene': 'AOC1', 'transcript': 'NM_001091.4', 'hgvsc': 'c.995C>T', 'hgvsp': 'p.Ser332Phe', 'build': 'hg38', 'location': 'chr7:150857465', 'chrom': 'chr7', 'pos': 150857465, 'ref': 'C', 'alt': 'T',
        'genotypes': {
            'C/C': {'zygosity': 'Homozigoto Referência', 'result': 'padrão', 
                    'interpretation': "O genótipo C/C (Ser/Ser) é o de referência (wild-type) e está associado à atividade enzimática normal da DAO."},
            'C/T': {'zygosity': 'Heterozigoto', 'result': 'minimamente reduzida ou normal', 
                    'interpretation': "O genótipo C/T (Ser/Phe) tem um efeito mínimo ou negligenciável na atividade da DAO."},
            'T/T': {'zygosity': 'Homozigoto Variante', 'result': 'minimamente reduzida', 
                    'interpretation': "O genótipo T/T (Phe/Phe) está associado a uma leve redução na atividade da DAO."}
        }
    },
    'rs1049793': {
        'gene': 'AOC1', 'transcript': 'NM_001091.4', 'hgvsc': 'c.1933C>G', 'hgvsp': 'p.His645Asp', 'build': 'hg38', 'location': 'chr7:150860577', 'chrom': 'chr7', 'pos': 150860577, 'ref': 'C', 'alt': 'G',
        'genotypes': {
            'C/C': {'zygosity': 'Homozigoto Referência', 'result': 'padrão', 
                    'interpretation': "O genótipo C/C (His/His) é o de referência (wild-type) e está associado à atividade enzimática normal da DAO."},
            'C/G': {'zygosity': 'Heterozigoto', 'result': 'reduzida', 
                    'interpretation': "O genótipo C/G (His/Asp) causa uma perda de atividade da DAO de aproximadamente ~34% (retendo ~66% da função normal)."},
            'G/G': {'zygosity': 'Homozigoto Variante', 'result': 'reduzida', 
                    'interpretation': "O genótipo G/G (Asp/Asp) causa uma perda de atividade da DAO de aproximadamente ~49% (retendo ~51% da função normal)."}
        }
    }
}

def genotype_position(bam_file, ref_file, chrom, pos, min_depth=30, min_base_quality=20, min_mapping_quality=20):
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
    """Gera o laudo HTML com formato de texto simples para compatibilidade com Word/RTF."""
    first_variant_info = next(iter(VARIANTS.values()))
    genome_build = first_variant_info['build']
    gene_name = first_variant_info['gene']
    transcript_id = first_variant_info['transcript']

    html_style = """
    <style>
        body { 
            font-family: "Times New Roman", Times, serif; 
            line-height: 1.6; 
            color: #333; 
            max-width: 900px; 
            margin: 40px auto; 
            font-size: 12px;
        }
        h1 { 
            font-size: 22px; 
            text-align: center; 
            color: #2c3e50; 
            border-bottom: 2px solid #3498db; 
            padding-bottom: 10px;
        }
        hr { 
            border: 0; 
            border-top: 1px solid #ddd; 
            margin: 40px 0 25px 0;
        }
        strong { color: #000; }
        .note { font-style: italic; }
    </style>
    """

    # --- LÓGICA DE PRÉ-PROCESSAMENTO PARA O RESUMO DINÂMICO ---
    genotype_results_for_interpretation = {}
    found_risk_variants = []
    max_severity = 0
    severity_map = {
        'padrão': 0, 'minimamente reduzida ou normal': 1, 'levemente reduzida': 2,
        'moderadamente reduzida': 3, 'reduzida': 4, 'reduzida (↓ ~34%)': 4,
        'severamente reduzida (↓ ~49%)': 5
    }
    for rsid, data in results.items():
        variant_info = VARIANTS[rsid]
        genotype_call = data['genotype']
        normalized_genotype = "/".join(sorted(genotype_call.split('/'))) if '/' in genotype_call else genotype_call
        genotype_results_for_interpretation[rsid] = {'call': genotype_call, 'key': None}
        key = ""
        if "Baixa Cobertura" not in normalized_genotype:
            ref, alt = variant_info['ref'], variant_info['alt']
            if normalized_genotype == f"{ref}/{ref}": key = f"{ref}/{ref}"
            elif normalized_genotype == f"{ref}/{alt}": key = f"{ref}/{alt}"
            elif normalized_genotype == f"{alt}/{alt}": key = f"{alt}/{alt}"
        if key:
            genotype_results_for_interpretation[rsid]['key'] = key
            info = variant_info['genotypes'][key]
            current_severity = severity_map.get(info['result'], 0)
            if current_severity > 0:
                variant_name_str = f"{variant_info['hgvsc']} ({variant_info['hgvsp']})"
                genotype_str = f"<strong>{key}</strong> ({info['zygosity']})"
                found_risk_variants.append({'name': variant_name_str, 'genotype': genotype_str})
                if current_severity > max_severity:
                    max_severity = current_severity

    # --- Construção do HTML ---
    html_content = f"""
    <!DOCTYPE html><html lang="pt-br"><head><meta charset="UTF-8"><title>Laudo Genético - {sample_id}</title>{html_style}</head><body>
    <h1>Deficiência de diamina oxidase (DAO)</h1>
    
    <hr style="margin-top: 25px;">
    <p>
        <strong>DADOS DO PACIENTE</strong><br><br>
        <strong>PACIENTE:</strong> {sample_id}<br>
        <strong>DATA DO LAUDO:</strong> {datetime.date.today().strftime('%d/%m/%Y')}<br>
        <strong>EXAME:</strong> Análise de variantes no gene AOC1 (DAO)<br>
        <strong>ARQUIVO ANALISADO:</strong> {os.path.basename(bam_filename)}<br>
    </p>

    <hr>
    <p>
        <strong>INFORMAÇÕES DA ANÁLISE</strong><br><br>
        <strong>Gene:</strong> {gene_name}<br>
        <strong>Transcrito:</strong> {transcript_id}<br>
        <strong>Referência genômica:</strong> {genome_build.upper()}<br>
        <br>
        <strong>Variantes Analisadas:</strong><br>
    """
    for rsid, v_info in VARIANTS.items():
        variant_line = f"{v_info['hgvsc']} ({v_info['hgvsp']})  |  ID: <strong>{rsid}</strong>  |  Posição: <strong>{v_info['location']}</strong>"
        html_content += f"{variant_line}<br>"
    html_content += "<br></p>"

    html_content += f'<hr><p><strong>RESULTADO</strong><br><br>'
    
    for rsid, data in results.items():
        variant_info = VARIANTS[rsid]
        genotype_call = data['genotype']
        variant_hgv = f"<strong>{variant_info['hgvsc']} ({variant_info['hgvsp']})</strong>"
        genotype_zygosity, activity, interpretation = "Indeterminado", "Inconclusivo", ""
        
        if "Baixa Cobertura" in genotype_call:
            genotype_zygosity = f"<strong>{genotype_call}</strong>"
            interpretation = "A análise foi <strong>inconclusiva</strong> devido à baixa cobertura de sequenciamento nesta posição, o que impede uma genotipagem confiável."
        else:
            key = genotype_results_for_interpretation[rsid]['key']
            if key:
                info = variant_info['genotypes'][key]
                genotype_zygosity = f"<strong>{key}</strong> ({info['zygosity']})"
                activity = f"<strong>{info['result']}</strong>"
                interpretation = info['interpretation']
            else:
                genotype_zygosity = f"<strong>{genotype_call}</strong> (Não canônico)"
                interpretation = "O resultado para esta variante foi inconclusivo. O genótipo observado não corresponde aos alelos de referência/alternativo esperados."
        
        result_line = f"{variant_hgv}  |  Genótipo: {genotype_zygosity}  |  Atividade Enzimática: {activity}"
        
        html_content += f"""
            {result_line}<br>
            <span style="font-style: italic;">{interpretation}</span><br><br>
        """
    html_content = html_content.rstrip('<br><br>') + '</p>'
    
    # --- Geração do resumo dinâmico ---
    summary_html = ""
    num_risk_variants = len(found_risk_variants)
    if num_risk_variants == 0:
        summary_html = "Nenhuma das variantes de risco analisadas foi detectada. O perfil genético do paciente é compatível com uma <strong>atividade normal</strong> da enzima DAO."
    else:
        formatted_variants_list = [f"<strong>{v['name']}</strong> com o genótipo {v['genotype']}" for v in found_risk_variants]
        if num_risk_variants == 1:
            verb, variant_text = "Foi detectada", "a variante"
            variant_list_str = formatted_variants_list[0]
        else:
            verb, variant_text = "Foram detectadas", "as variantes"
            variant_list_str = ", ".join(formatted_variants_list[:-1]) + " e " + formatted_variants_list[-1]
        
        if max_severity <= 2: conclusion_text = "sugerem que o impacto na atividade enzimática é <strong>mínimo ou leve</strong>."
        elif max_severity == 3: conclusion_text = "sugerem uma <strong>atividade enzimática moderadamente reduzida</strong>."
        elif max_severity == 4: conclusion_text = "sugerem uma <strong>atividade enzimática reduzida</strong>."
        else: conclusion_text = "sugerem uma <strong>atividade enzimática severamente reduzida</strong>."
        
        summary_html = f"{verb} {variant_text} {variant_list_str}, que, em conjunto, {conclusion_text}"
    
    # --- Bloco de texto fixo ---
    fixed_interpretation_text = """
    <br><br>
    Portadores das variantes p.(Thr16Met), p.(Ser332Phe) e p.(His645Asp) mostram atividade enzimática
    diminuída em comparação aos não portadores das referidas variantes (1,2). Foi descrito que a variante
    p.(Thr16Met) tem efeito intermediário e a variante p.(His645Asp) tem efeito mais severo, sendo a
    variante p.(Ser332Phe) a que afeta menos a atividade enzimática.<br><br>
    A variante c.47C>T p.(Thr16Met) foi associada a um maior risco de hipersensibilidade a anti-
    inflamatórios não esteroides (3).<br><br>
    As variantes de suscetibilidade são aquelas alterações no DNA que afetam a atividade de uma
    enzima/metabólito/proteína influenciando o risco de desenvolver uma doença. Para que esta variante
    genética se expresse muitas vezes, requer o envolvimento de fatores não genéticos.<br><br>
    A deficiência de DAO é uma alteração no metabolismo da histamina alimentar devido a uma menor
    atividade da enzima diamina oxidase, produzindo um acúmulo de histamina no plasma. O quadro clínico
    varia desde assintomáticos até pacientes com enxaquecas, distúrbios gastrointestinais, dermatológicos,
    entre outros.<br><br>
    <strong>Este laudo deve ser interpretado por um especialista dentro do contexto clínico e da história familiar do
    paciente, juntamente com outros achados laboratoriais. Caso o médico solicitante julgue necessário,
    recomenda-se aconselhamento genético.</strong>
    """
    
    # --- Montagem da seção de Interpretação final ---
    html_content += f'<hr><p><strong>INTERPRETAÇÃO DO RESULTADO</strong><br><br>{summary_html}{fixed_interpretation_text}</p>'

    html_content += f'<hr><p><strong>CONTROLE DE QUALIDADE DA GENOTIPAGEM</strong><br><br>'
    for rsid, data in results.items():
        variant_info = VARIANTS[rsid]
        coverage_str = f"{data['coverage']}x"
        counts_str = ", ".join([f"{base} {count}x" for base, count in data['counts'].items()]) or "N/A"
        quality_line = f"{variant_info['hgvsc']} ({variant_info['hgvsp']})  |  Cobertura (BQ≥{min_bq}): <strong>{coverage_str}</strong>  |  Leituras: <strong>{counts_str}</strong>"
        html_content += f"{quality_line}<br>"
    html_content += '<br></p>'

    html_content += f"""
        <hr>
        <p>
            <strong>TÉCNICA APLICADA</strong><br><br>
            A análise foi realizada a partir de DNA genômico extraído de forma automatizada (Maxwell®, Promega). O preparo da biblioteca de DNA foi conduzido com a tecnologia de captura em alvo "xGen™ DNA Library Prep EZ" e o painel de enriquecimento "xGen™ Exome Hyb Panel v2" (Integrated DNA Technologies - IDT), que abrange as regiões codificantes e intrônicas adjacentes (±10 pares de base) dos genes nucleares. O sequenciamento de nova geração (NGS) foi executado em uma plataforma MGI Tech DNABSEQ G400.<br><br>
            O processamento primário dos dados brutos, incluindo o alinhamento das leituras ao genoma humano de referência (GRCh38/hg38) e a chamada de variantes, foi realizado através do pipeline Dragen Enrichment (v. 4.2.4).<br><br>
            As variantes de interesse no gene AOC1 foram submetidas a uma análise de genotipagem focada, executada por um pipeline próprio que consulta diretamente o arquivo de alinhamento (BAM) para verificar a cobertura e a identidade das bases em cada posição genômica alvo.<br><br>
            A interpretação dos genótipos é baseada na base de conhecimento interno compilado a partir das referências científicas listadas abaixo.<br><br>
            <strong>Limitações e Parâmetros de Qualidade:</strong> A análise se restringe apenas às variantes descritas neste laudo. A análise não determina a fase haplotípica (configuração <i>cis/trans</i>) entre as variantes identificadas. A genotipagem depende da cobertura efetiva na posição de interesse (limite mínimo: 30x).<br><br>
            <span class="note"><strong>Nota sobre Qualidade (BQ):</strong> A genotipagem realizada considera apenas leituras com Qualidade de Base (Base Quality - BQ) igual ou superior a 20.<br>A BQ é uma medida de confiança na identificação correta de cada base nucleotídica (A, C, G, T) pela plataforma de sequenciamento.</span><br>
        </p>

        <hr>
        <p>
            <strong>REFERÊNCIAS BIBLIOGRÁFICAS</strong><br><br>
            1) Maintz L, et al. Association of single nucleotide polymorphisms in the diamine oxidase gene with diamine oxidase serum activities. Allergy. 2011 Jul;66(7):893-902.<br>
            2) Ayuso P, et al. Genetic variability of human diamine oxidase: occurrence of three nonsynonymous polymorphisms and study of their effect on serum enzyme activity. Pharmacogenet Genomics. 2007 Sep;17(9):687-93.<br>
            3) Agúndez José A G et al. The diamine oxidase gene is associated with hypersensitivity response to non-steroidal anti-inflammatory drugs. PLoS One. 2012;7(11):e47571.<br>
        </p>
    </body></html>
    """
    return html_content

def main():
    parser = argparse.ArgumentParser(description="Gera um laudo para variantes do gene AOC1 a partir de um arquivo BAM.")
    parser.add_argument("bam_file", help="Caminho para o arquivo BAM do paciente.")
    parser.add_argument("ref_file", help="Caminho para o arquivo FASTA do genoma de referência.")
    parser.add_argument("--sample_id", help="ID do paciente/amostra para o laudo.", default="Amostra Anônima")
    parser.add_argument("--output_file", help="Nome do arquivo HTML para salvar o laudo (ex: laudo.html).", default="laudo_aoc1.html")
    parser.add_argument("--min_depth", type=int, default=30, help="Profundidade mínima de cobertura para genotipagem. Padrão: 30x.")
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
            min_depth=args.min_depth,
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