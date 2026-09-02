#!/usr/bin/env python3
"""
generate_final_shortlists.py — Final Ranked Shortlists for Track 3 Submission
==============================================================================
Generates 5 submission outputs using n=20 WGCNA data with explicit ranking functions.

OUTPUTS (written to results/shortlist/):
  1. host_genes.csv           — ≤10 ranked host genes
  2. microbial_taxa.csv        — 5 ranked microbial taxa
  3. pathways.csv              — 3 ranked host pathways
  4. network_edges.csv         — Final edge list for integration figure
  5. ranking_methodology.md    — Full ranking documentation

RANKING FORMULAS:
  Host genes:  Score = |partial_r(WG|sex)| × |kME| × growth_module_bonus
  Microbial:   Score = |CLR_fold_diff| × (1−contamination_score) × direction_multiplier
  Pathways:    Score = edge_count × |module_partial_r| × biological_relevance_factor

CRITICAL: All microbial associations are directional hypotheses only (n=1 pooled).
No cross-layer covariance is claimed. Statistics from n=20 WGCNA.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[3]


def load_data():
    """Load all input data, including v3 stability reports."""
    print("LOADING DATA...")
    me = pd.read_csv(BASE / 'data/processed/wgcna/me_matrix.tsv', sep='\t', index_col=0)
    hub = pd.read_csv(BASE / 'data/processed/wgcna/hub_genes.tsv', sep='\t')
    modules = pd.read_csv(BASE / 'results/tables/wgcna_modules.csv', sep='\t')
    meta = pd.read_csv(BASE / 'data/raw/sra/PRJNA875278/metadata.tsv', sep='\t', index_col=0)
    taxa = pd.read_csv(BASE / 'data/processed/clr_profiles/taxa_direction.tsv', sep='\t')
    edges = pd.read_csv(BASE / 'results/shortlist/network_edges.csv')
    
    # v3: Load stability reports if they exist
    taxon_conf = None
    wgcna_stability = None
    deseq2_modules = None
    conf_path = BASE / 'results/reports/taxon_confidence_report.tsv'
    wgcna_path = BASE / 'results/reports/wgcna_module_stability_summary.tsv'
    deseq2_path = BASE / 'results/reports/deseq2_module_eigengene_test.tsv'
    if conf_path.exists():
        taxon_conf = pd.read_csv(conf_path, sep='\t')
        print(f"  Taxon confidence report: {len(taxon_conf)} taxa")
    if wgcna_path.exists():
        wgcna_stability = pd.read_csv(wgcna_path, sep='\t')
        print(f"  WGCNA module stability: {len(wgcna_stability)} modules")
    if deseq2_path.exists():
        deseq2_modules = pd.read_csv(deseq2_path, sep='\t')
        print(f"  DESeq2 module eigengene test: {len(deseq2_modules)} modules")
    
    print(f"  WGCNA: {len(me)} samples, {len(hub)} hub genes, {len(modules)} gene assignments")
    print(f"  Taxa: {len(taxa)}  |  Edges: {len(edges)}  |  Metadata: {len(meta)} samples")
    return me, hub, modules, meta, taxa, edges, taxon_conf, wgcna_stability, deseq2_modules


def module_stats(me, meta):
    """
    Compute module-trait partial correlations (n=20).
    Returns DataFrame sorted by |partial_r(WG|sex)|.
    """
    common = sorted(set(me.index) & set(meta.index))
    wg = meta.loc[common, 'weight_gain'].values.astype(float)
    sex = np.array([1.0 if str(meta.loc[s, 'sex']).lower() == 'male' else 0.0
                    for s in common])
    r_wg_sex = np.corrcoef(wg, sex)[0, 1]

    rows = []
    for col in sorted(c for c in me.columns if c.startswith('ME')):
        mn = col.replace('ME', '')
        mv = me.loc[common, col].values.astype(float)
        rw = np.corrcoef(mv, wg)[0, 1]
        rs = np.corrcoef(mv, sex)[0, 1]
        num = rw - rs * r_wg_sex
        denom = np.sqrt((1 - rs**2) * (1 - r_wg_sex**2))
        rp = num / denom if denom != 0 else 0.0
        rows.append({
            'module': f'M{mn}',
            'r_wg': round(rw, 4),
            'r_sex': round(rs, 4),
            'partial_r': round(rp, 4),
            'abs_partial_r': round(abs(rp), 4),
        })
    return pd.DataFrame(rows).sort_values('abs_partial_r', ascending=False)


def rank_host_genes(hub, me, meta, modules, wgcna_stability=None, deseq2_modules=None, max_n=10):
    """
    Score = |partial_r(WG|sex)| × |kME| × growth_module_bonus × module_stability_bonus
            × deseq2_concordance_bonus (v3.1).
    Max 3 genes per module for diversity.

    v3.1: DESeq2 concordance bonus from module-level eigengene LM test.
    ✓ concordant → 1.0, ✗ discordant → 0.5 (penalized).
    M7 specifically penalized if DESeq2 contradicts.
    """
    print("\n=== RANKING HOST GENES (target: <=10) ===")
    mstats = module_stats(me, meta)
    top5 = set(mstats.head(5)['module'])
    print(f"  Top 5 growth modules (by |partial r|): {sorted(top5)}")
    for _, r in mstats.head(5).iterrows():
        ng = len(modules[modules['module'] == r['module']])
        print(f"    {r['module']}: partial_r={r['partial_r']:+.4f}, "
              f"r_wg={r['r_wg']:+.4f}, r_sex={r['r_sex']:+.4f}, {ng} genes")

    merged = hub.merge(
        mstats[['module', 'partial_r', 'abs_partial_r']], on='module', how='left')
    merged['bonus'] = merged['module'].isin(top5).map({True: 1.0, False: 0.5})

    # v3: Module stability bonus (from bootstrap)
    stability_bonus = {}
    if wgcna_stability is not None:
        for _, row in wgcna_stability.iterrows():
            mod = row['module']
            tier = row.get('module_stability_tier', 'MEDIUM')
            bonus = {'HIGH': 1.0, 'MEDIUM': 0.8, 'LOW': 0.5}.get(tier, 0.7)
            stability_bonus[mod] = bonus
        print(f"  v3: Module stability bonuses loaded for {len(stability_bonus)} modules")
    merged['stability_bonus'] = merged['module'].map(
        lambda m: stability_bonus.get(m, 0.7))

    # v3.1: DESeq2 module-level concordance bonus/penalty
    deseq2_concordance = {}
    if deseq2_modules is not None:
        for _, row in deseq2_modules.iterrows():
            mod = row['module']
            is_conc = row.get('direction_concordant', True)
            pval = row.get('wg_pvalue', 0.5)
            coef = row.get('wg_coef', 0.0)
            # DESeq2 support for WGCNA finding:
            # AGREE_SIGNIFICANT: concordant + p<0.1 → boost 1.2
            # AGREE: concordant + 0.1<=p<=0.5 → neutral 1.0
            # WEAK_EFFECT: concordant but p>0.5 (negligible WG effect) → penalize 0.5
            # DISAGREE: discordant → penalize 0.5
            # STRONG_DISAGREE: discordant + p<0.1 → heavy penalty 0.3
            if not is_conc and pval < 0.1:
                deseq2_concordance[mod] = (0.3, 'STRONG_DISAGREE')
            elif not is_conc:
                deseq2_concordance[mod] = (0.5, 'DISAGREE')
            elif is_conc and pval < 0.1:
                deseq2_concordance[mod] = (1.2, 'AGREE_SIGNIFICANT')
            elif is_conc and pval > 0.5:
                deseq2_concordance[mod] = (0.5, 'WEAK_EFFECT')
            else:
                deseq2_concordance[mod] = (1.0, 'AGREE')
        
        n_agree_sig = sum(1 for v in deseq2_concordance.values() if v[1] == 'AGREE_SIGNIFICANT')
        n_weak = sum(1 for v in deseq2_concordance.values() if v[1] == 'WEAK_EFFECT')
        n_disagree = sum(1 for v in deseq2_concordance.values() if v[1] in ('DISAGREE', 'STRONG_DISAGREE'))
        print(f"  v3.1: DESeq2 module concordance: {n_agree_sig} sig-agree, "
              f"{n_weak} weak-effect, {n_disagree} disagree")
    merged['deseq2_bonus'] = merged['module'].map(
        lambda m: deseq2_concordance.get(m, (1.0, 'UNKNOWN'))[0])
    merged['deseq2_status'] = merged['module'].map(
        lambda m: deseq2_concordance.get(m, (1.0, 'UNKNOWN'))[1])

    merged['score'] = (
        merged['abs_partial_r'].fillna(0) *
        merged['kME'].abs() *
        merged['bonus'] *
        merged['stability_bonus'] *
        merged['deseq2_bonus']
    )
    merged = merged.sort_values('score', ascending=False)

    picked, counts = [], {}
    for _, r in merged.iterrows():
        m = r['module']
        if counts.get(m, 0) >= 3:
            continue
        if len(picked) >= max_n:
            break
        picked.append(r)
        counts[m] = counts.get(m, 0) + 1

    result = pd.DataFrame(picked)
    result.insert(0, 'rank', range(1, len(result) + 1))
    for _, r in result.iterrows():
        print(f"  #{int(r['rank']):2d} {r['gene']:<25s} M{r['module']:<5s} "
              f"kME={r['kME']:+.4f} partial_r={r['partial_r']:+.4f} "
              f"score={r['score']:.6f} deseq2={r.get('deseq2_status', 'N/A')}")
    return result, mstats


def rank_microbial_taxa(taxa, taxon_conf=None, max_n=5):
    """
    Score = |CLR_fold_diff| × (1-contamination_score) × direction_multiplier
            × stability_multiplier (v3).
    Exclude HIGH contamination risk and 'Unknown' taxon.
    
    v3: Incorporates combined_stability from pseudocount+LOTO analysis.
    HIGH stability → 1.2 (boost), MEDIUM → 1.0, LOW → 0.7 (penalize).
    """
    print("\n=== RANKING MICROBIAL TAXA (target: 5) ===")
    df = taxa.copy()

    if 'abs_fold_diff' not in df.columns:
        df['abs_fold_diff'] = df['fold_diff'].abs()
    if 'contamination_score' not in df.columns:
        df['contamination_score'] = df['contaminant_flag'].astype(float)

    # v3: Merge stability info if available
    stability_mult = {}
    if taxon_conf is not None:
        for _, row in taxon_conf.iterrows():
            tier = row.get('stability_tier', 'MEDIUM')
            mult = {'HIGH': 1.2, 'MEDIUM': 1.0, 'LOW': 0.7}.get(tier, 1.0)
            stability_mult[row['taxon']] = (mult, tier)
        print(f"  v3: Taxon stability data loaded for {len(stability_mult)} taxa")
    
    df['stability_mult'] = df['taxon'].map(
        lambda t: stability_mult.get(t, (1.0, 'UNKNOWN'))[0])
    df['stability_tier'] = df['taxon'].map(
        lambda t: stability_mult.get(t, (1.0, 'UNKNOWN'))[1])

    df['dir_mult'] = df['direction'].apply(
        lambda d: 1.0 if 'Jumper' in str(d) else (0.7 if 'Laggard' in str(d) else 0.0))
    df['score'] = (
        df['abs_fold_diff'].fillna(0) *
        (1 - df['contamination_score'].fillna(0)) *
        df['dir_mult'] *
        df['stability_mult'].fillna(1.0)
    )

    mask = ((df.get('contamination_risk', 'LOW') != 'HIGH') &
            (df['taxon'].str.lower() != 'unknown'))
    df_filt = df[mask].sort_values('score', ascending=False)

    picked = df_filt.head(max_n).copy()
    picked.insert(0, 'rank', range(1, len(picked) + 1))
    for _, r in picked.iterrows():
        print(f"  #{int(r['rank'])} {r['taxon']:<40s} "
              f"CLR_fd={r['fold_diff']:+.4f} dir={r['direction']:<20s} "
              f"score={r['score']:.4f} contam_score={r['contamination_score']:.2f} "
              f"stability={r.get('stability_tier', 'N/A')}")
    return picked


def rank_pathways(edges, me, meta, max_n=3):
    """
    Score = edge_count × biological_relevance_factor × |top_module_partial_r|.
    """
    print("\n=== RANKING PATHWAYS (target: 3) ===")
    func_e = edges[edges['edge_basis'] == 'predicted_function_overlap']
    pw_counts = func_e.groupby('target_name').size().reset_index(name='edge_count')

    relevance = {
        'Amino sugar and nucleotide sugar metabolism': 1.5,
        'Lysosome': 1.5,
        'Citrate cycle (TCA cycle)': 1.0,
        'Pyruvate metabolism': 1.0,
        'Glyoxylate and dicarboxylate metabolism': 0.9,
        'PPAR signaling pathway': 1.0,
        'Fatty acid metabolism': 0.9,
        'Glycerolipid metabolism': 0.7,
    }
    pw_counts['bio_rel'] = pw_counts['target_name'].map(
        lambda x: relevance.get(x, 0.7))

    mstats = module_stats(me, meta)
    top_r = mstats.iloc[0]['abs_partial_r']
    pw_counts['score'] = pw_counts['edge_count'] * pw_counts['bio_rel'] * top_r

    picked = pw_counts.sort_values('score', ascending=False).head(max_n).copy()
    picked.insert(0, 'rank', range(1, len(picked) + 1))

    for i, row in picked.iterrows():
        sel = func_e[func_e['target_name'] == row['target_name']]
        if 'source_taxa' in sel.columns:
            taxa = set()
            for v in sel['source_taxa'].dropna():
                for part in str(v).split('; '):
                    taxa.add(part.split(' (')[0])
            picked.at[i, 'source_taxa'] = '; '.join(sorted(t for t in taxa if t))
        else:
            picked.at[i, 'source_taxa'] = '; '.join(sorted(sel['source_name'].unique()))

    for _, r in picked.iterrows():
        print(f"  #{int(r['rank'])} {r['target_name']:<55s} "
              f"edges={r['edge_count']:3d} bio_rel={r['bio_rel']:.1f} "
              f"score={r['score']:.4f}")
    return picked, mstats


def write_outputs(host, taxa_r, pathways, mstats, edges, modules, out_dir):
    """Write all 5 submission output files."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1. Host genes ----
    hrows = []
    for _, g in host.iterrows():
        ms = mstats[mstats['module'] == g['module']].iloc[0]
        ng = len(modules[modules['module'] == g['module']])
        hrows.append({
            'rank': int(g['rank']),
            'gene_id': g['gene'],
            'source_analysis': 'WGCNA co-expression network (n=20, power=15, R^2=0.939)',
            'statistic_attributable_to': (
                f"kME={g['kME']:.4f} (Pearson r of gene with module eigengene, n=20). "
                f"Module {g['module']} partial r(WG|sex)={ms['partial_r']:+.4f} (n=20). "
                f"Module size: {ng} genes."
            ),
            'associated_module': g['module'],
            'kME': round(g['kME'], 4),
            'module_partial_r_wg_given_sex': ms['partial_r'],
            'module_r_wg': ms['r_wg'],
            'module_r_sex': ms['r_sex'],
            'trait': 'weight_gain',
            'n_rnaseq_libraries': 20,
            'sex_confound_note': (
                'Weight gain and sex confounded (r=-0.762). '
                'Partial correlation controls for sex. '
                'Pre-WGCNA PC removal (Fix 3) removed 66.2% confounded variance.'
            ),
            'rank_score': round(g['score'], 6),
            'ranking_formula': '|partial_r(WG|sex)| x |kME| x growth_module_bonus x module_stability_bonus (v3)',
            'module_stability_bonus': g.get('stability_bonus', 'N/A'),
            'evidence_basis': 'WGCNA module-trait partial correlation (n=20)',
        })
    pd.DataFrame(hrows).to_csv(out_dir / 'host_genes.csv', index=False)
    print(f"\n  [OK] host_genes.csv: {len(hrows)} genes")

    # ---- 2. Microbial taxa ----
    mrows = []
    for _, r in taxa_r.iterrows():
        mrows.append({
            'rank': int(r['rank']),
            'row_label_in_supplied_table': r['taxon'],
            'approach_used': 2,
            'statistic_attributable_to': (
                f"CLR fold-diff = {r['fold_diff']:.4f}. "
                f"Raw reads: Jumper={int(r['Jumper'])}, Laggard={int(r['Laggard'])}. "
                f"Contamination score: {r['contamination_score']:.2f} (0=clean). "
                f"Stability tier: {r.get('stability_tier', 'N/A')}. "
                f"NO statistical test -- n=1 pooled per group."
            ),
            'direction_of_effect': r['direction'],
            'jumper_raw_reads': int(r['Jumper']),
            'laggard_raw_reads': int(r['Laggard']),
            'clr_fold_diff': round(r['fold_diff'], 4),
            'contaminant_risk': r.get('contamination_risk', 'LOW'),
            'contamination_score': r.get('contamination_score', 0.0),
            'contaminant_rationale': r.get('contamination_rationale', ''),
            'stability_tier': r.get('stability_tier', 'N/A'),
            'stability_multiplier': r.get('stability_mult', 'N/A'),
            'mechanism_summary': (
                f"{r['taxon']} is {r['direction']} in Jumper vs Laggard pooled 16S. "
                f"CLR fold-diff: {r['fold_diff']:.2f}. "
                f"Cannot be statistically tested (n=1 pooled). "
                f"Contamination risk: {r.get('contamination_risk', 'LOW')}. "
                f"Stability tier: {r.get('stability_tier', 'N/A')}."
            ),
            'evidence_basis': 'CLR direction from pooled ASV (n=1/group) -- directional hypothesis only',
            'n_microbial_pools': 2,
            'rank_score': round(r['score'], 6),
            'ranking_formula': '|CLR_fold_diff| x (1-contamination_score) x direction_multiplier x stability_multiplier (v3)',
        })
    pd.DataFrame(mrows).to_csv(out_dir / 'microbial_taxa.csv', index=False)
    print(f"  [OK] microbial_taxa.csv: {len(mrows)} taxa")

    # ---- 3. Pathways ----
    prows = []
    for _, r in pathways.iterrows():
        bio_note = (
            'Chitin degradation -> GlcNAc -> exoskeleton synthesis (crustacean molting)'
            if r['bio_rel'] >= 1.5 else
            'SCFA (butyrate/propionate) sensing -> host energy metabolism'
            if r['target_name'] == 'PPAR signaling pathway' else
            'Central energy metabolism -- organic acid fermentation -> host TCA cycle'
            if r['bio_rel'] >= 1.0 else
            'Ancillary pathway with plausible microbial modulation'
        )
        prows.append({
            'rank': int(r['rank']),
            'target_name': r['target_name'],
            'n_microbial_sources': r['edge_count'],
            'source_taxa': r['source_taxa'],
            'edge_basis_types': 'predicted_function_overlap (KEGG GENOME-derived pathway index + curated host bridge)',
            'top_wgcna_module': mstats.iloc[0]['module'],
            'top_module_partial_r_wg': mstats.iloc[0]['partial_r'],
            'biological_relevance_factor': r['bio_rel'],
            'biological_rationale': bio_note,
            'rank_score': round(r['score'], 6),
            'ranking_formula': 'edge_count x |module_partial_r| x biological_relevance_factor',
            'evidence_basis': 'KEGG GENOME-derived microbial pathway index -> curated host pathway bridge (no cross-layer correlation; n=1 pools descriptive)',
        })
    pd.DataFrame(prows).to_csv(out_dir / 'pathways.csv', index=False)
    print(f"  [OK] pathways.csv: {len(prows)} pathways")

    # ---- 4. Network edges (final figure) ----
    ranked_taxa_names = set(taxa_r['taxon'].values)
    ranked_pw_names = set(pathways['target_name'].values)
    final_e = edges[
        edges['source_name'].isin(ranked_taxa_names) |
        edges['target_name'].isin(ranked_pw_names)
    ].copy()

    final_e['edge_meaning'] = final_e['edge_basis'].map({
        'predicted_function_overlap':
            'KEGG GENOME-derived microbial pathway index -> curated host KEGG '
            'pathway bridge. NO cross-layer correlation claimed. '
            'Direction descriptive only (n=1 pooled).',
        'phenotype_concordance':
            'Taxon enriched in Jumper pooled sample vs Laggard. '
            'Directional hypothesis only. NOT a correlation.'
    })
    final_e['n_rnaseq_libraries_supporting'] = 20
    final_e['n_microbial_pools'] = 2
    final_e.to_csv(out_dir / 'network_edges.csv', index=False)
    print(f"  [OK] network_edges.csv: {len(final_e)} edges")

    # ---- 5. Ranking methodology ----
    top_mods_str = "\n".join(
        f"  - **{r['module']}**: partial r(WG|sex)={r['partial_r']:+.4f}, "
        f"r(WG)={r['r_wg']:+.4f}, r(sex)={r['r_sex']:+.4f}"
        for _, r in mstats.head(5).iterrows()
    )
    host_str = "\n".join(
        f"  - **#{int(r['rank'])} {r['gene']}** (M{r['module']}): "
        f"kME={r['kME']:+.4f}, score={r['score']:.6f}, "
        f"DESeq2={r.get('deseq2_status', 'N/A')}"
        for _, r in host.iterrows()
    )
    taxa_str = "\n".join(
        f"  - **#{int(r['rank'])} {r['taxon']}**: "
        f"CLR fold-diff={r['fold_diff']:+.4f}, {r['direction']}, "
        f"contam_score={r['contamination_score']:.3f}, "
        f"stability={r.get('stability_tier', 'N/A')}, "
        f"score={r['score']:.6f}"
        for _, r in taxa_r.iterrows()
    )
    pw_str = "\n".join(
        f"  - **#{int(r['rank'])} {r['target_name']}**: "
        f"{r['edge_count']} microbial edges, bio_relevance={r['bio_rel']:.1f}, "
        f"score={r['score']:.6f}"
        for _, r in pathways.iterrows()
    )

    with open(out_dir / 'ranking_methodology.md', 'w') as f:
        f.write(f"""# Ranking Methodology -- Track 3 Host-Microbe Integration (v3.1, n=20)

## 1. Host Gene Ranking

**Formula**: `Score = |partial_r(WG|sex)| × |kME| × growth_module_bonus × module_stability_bonus × deseq2_concordance_bonus`

| Component | Source | n | v3.1 Bonus |
|-----------|--------|---|------------|
| partial_r(WG\\|sex) | Partial correlation of module eigengene with weight_gain, controlling for sex | 20 RNA-seq libraries | — |
| kME | Module membership: Pearson r(gene, module eigengene) | 20 samples | — |
| growth_module_bonus | 1.0 if module in top 5 by \\|partial_r\\|, else 0.5 | — | — |
| module_stability_bonus | Bootstrap stability tier (HIGH=1.0, MEDIUM=0.8, LOW=0.5) | 1000 bootstraps | v3 |
| deseq2_concordance_bonus | DESeq2 module eigengene LM test (AGREE_SIGNIFICANT=1.2, AGREE=1.0, WEAK_EFFECT=0.5, DISAGREE=0.5, STRONG_DISAGREE=0.3) | 20 samples, design=~sex+tissue+WG | v3.1 |

**Constraint**: Max 3 genes per module for biological diversity.

**Top 5 growth modules (n=20)**:
{top_mods_str}

**Selected host genes**:
{host_str}

---

## 2. Microbial Taxon Ranking

**Formula**: `Score = |CLR_fold_diff| × (1 − contamination_score) × direction_multiplier × stability_multiplier`

| Component | Source | n | v3.1 |
|-----------|--------|---|------|
| CLR_fold_diff | CLR(Jumper) − CLR(Laggard) | 2 pooled samples | — |
| contamination_score | Log-normal mixture model on read counts (Fix 1) | 191 taxa | — |
| direction_multiplier | 1.0 (Jumper-associated), 0.7 (Laggard-associated) | — | — |
| stability_multiplier | Pseudocount+LOTO combined stability tier (HIGH=1.2, MEDIUM=1.0, LOW=0.7) | 11 CLR variants × 191 LOTO iterations | v3 |

**Exclusions**: HIGH contamination risk (score ≥0.85), "Unknown" taxon.

**CRITICAL**: All microbial associations are DIRECTIONAL HYPOTHESES ONLY.
No p-values, FDR, or statistical tests. n=1 pooled per group.

**Selected microbial taxa**:
{taxa_str}

---

## 3. Pathway Ranking

**Formula**: `Score = edge_count × |module_partial_r| × biological_relevance_factor`

| Component | Source | n |
|-----------|--------|---|
| edge_count | Number of microbial KEGG pathways with function-overlap edges to this pathway | KEGG GENOME profiles |
| module_partial_r | \\|Partial r(WG\\|sex)\\| of strongest growth module | 20 samples |
| biological_relevance | 1.5 (chitin/amino-sugar), 1.0 (TCA/energy), 0.7 (other) | Expert curation |

**Selected pathways**:
{pw_str}

---

## 4. Validation Layers (v3.1)

### 4.1 Microbial Taxon Stability

| Analysis | Method | Key Finding |
|----------|--------|-------------|
| Pseudocount sensitivity | 11 CLR variants (pseudocount 0.1–50.0 + multiplicative + Bayesian) | 0/191 taxa pseudocount-stable; 179/191 flip direction |
| Leave-one-taxon-out (LOTO) | 191 jackknife iterations | Top 3 taxa LOTO-stable (1.000); C. koseri & K. variicola LOTO-unstable (0.000) |

### 4.2 Host Module Stability

| Analysis | Method | Key Finding |
|----------|--------|-------------|
| Bootstrap module-trait | 1000 resamples of n=20 with replacement | No module has 95% CI excluding zero; M6 top5_rate=0.748 |
| DESeq2 module eigengene LM | lm(ME ~ sex + tissue + WG), all 20 samples | M6: p=0.03 (significant); M7: p=0.81 (WG explains nothing) |
| DESeq2 per-gene cross-check | design=~sex+tissue+WG, continuous n=20 | 7/10 hub genes direction-concordant (kME-aware); 3/10 nominally significant |

### 4.3 Concordance Summary

| Layer | Concordance | Binomial p | Assessment |
|-------|-------------|------------|------------|
| Hub gene direction (v3.1 fixed) | 7/10 (70%) | p=0.17 | Not significant at α=0.05 with n=10 genes |
| Module eigengene direction | 13/20 (65%) | p=0.13 | Not significant at α=0.05 with n=20 modules |
| M6 (strongest module) | ✓ Concordant | LM p=0.03 | **Significant** — strongest validated finding |

---

## 5. What Each Edge Means

| Edge Type | Meaning |
|-----------|---------|
| `predicted_function_overlap` | KEGG GENOME-derived microbial pathway index → curated host pathway bridge. **NO cross-layer correlation.** |
| `phenotype_concordance` | Taxon more abundant in Jumper pool. **Directional hypothesis only.** |

## 6. Confound Handling

- **Sex-WG confound**: r(WG, sex) = −0.762 (females larger).
- **Post-hoc**: Partial correlation r(WG\\|sex) on module eigengenes.
- **Pre-WGCNA (Fix 3)**: PCA removal of sex/tissue PCs → 66.2% variance removed.
- **v3.1**: DESeq2 design includes sex + tissue as covariates; module LM confirms sex/tissue separation.
""")
    print(f"  [OK] ranking_methodology.md")

    return final_e


# =====================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("FINAL SUBMISSION SHORTLISTS -- Track 3 Host-Microbe Integration")
    print(f"v3: Stability-weighted ranking (pseudocount+LOTO+bootstrap+DESeq2)")
    print("=" * 70)

    me, hub, modules, meta, taxa, edges, taxon_conf, wgcna_stability, deseq2_modules = load_data()
    host_ranked, mstats = rank_host_genes(
        hub, me, meta, modules, wgcna_stability=wgcna_stability,
        deseq2_modules=deseq2_modules, max_n=10)
    taxa_ranked = rank_microbial_taxa(taxa, taxon_conf=taxon_conf, max_n=5)
    pathways_ranked, _ = rank_pathways(edges, me, meta, max_n=3)

    out_dir = BASE / 'results' / 'shortlist'
    final_edges = write_outputs(
        host_ranked, taxa_ranked, pathways_ranked, mstats,
        edges, modules, out_dir
    )

    print("\n" + "=" * 70)
    print("COMPLETE -- All files in results/shortlist/")
    print(f"  host_genes.csv:     {len(host_ranked)} genes (target <=10)")
    print(f"  microbial_taxa.csv: {len(taxa_ranked)} taxa (target 5)")
    print(f"  pathways.csv:       {len(pathways_ranked)} pathways (target 3)")
    print(f"  network_edges.csv:  {len(final_edges)} edges")
    print(f"  ranking_methodology.md")
    print("=" * 70)
