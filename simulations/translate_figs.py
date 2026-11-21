import re,sys
R={
'fig56.py':[
("ax[0].set_yticklabels(order)","ax[0].set_yticklabels(['Weekly ground, naive','Weekly ground, smart','Ground every 3 days, smart','Conventional autonomous','ΔF autonomous (relational)'])"),
("'mínimo físico (arrasto)'","'physical minimum (drag)'"),
("'Δv em regime (m/s por ano)'","'steady-state Δv (m/s per year)'"),
("'erro de altitude média, RMS (m)'","'mean-altitude error, RMS (m)'"),
("'ciclos de solo por ano'","'ground cycles per year'"),
("'Objetivo 1 — LEO 350 km, 60 dias × 3 sementes: autonomia atinge o mínimo físico com altitude 6–11× mais justa e zero ciclos de solo'","'Objective 1 — LEO 350 km, 60 days × 3 seeds: autonomy reaches the physical minimum with 6–11× tighter altitude and zero ground cycles'"),
("'Órbita a 1,75 km (período 58 h)'","'Orbit at 1.75 km (period 58 h)'"),
("'Sobrevoo hiperbólico a 7,25 km (tipo OSIRIS-REx)'","'Hyperbolic flyby at 7.25 km (OSIRIS-REx-like)'"),
("'Relacional M = ⟨h⟩²/(G·p_geom), SRP compensada'","'Relational M = ⟨h⟩²/(G·p_geom), SRP compensated'"),
("'OD em lote (μ + 6 estados, SRP modelada)'","'Batch OD (μ + 6 states, SRP modelled)'"),
("'Relacional sem compensar SRP'","'Relational, SRP not compensated'"),
("'1 % (requisito OSIRIS-REx p/ inserção)'","'1 % (OSIRIS-REx insertion requirement)'"),
("'duração do arco de dados'","'data-arc duration'"),
("'erro da massa, 1σ (%)'","'mass error, 1σ (%)'"),
("'Objetivo 2 — corpo tipo Bennu (GM 4,89 m³/s²), SRP 1×10⁻⁷ m/s², 40 Monte Carlo por ponto'","'Objective 2 — Bennu-like body (GM 4.89 m³/s²), SRP 1×10⁻⁷ m/s², 40 Monte Carlo runs per point'"),
],
'fig7.py':[
("'Relacional (SRP compensada)'","'Relational (SRP compensated)'"),("'OD massa pontual'","'Point-mass OD'"),("'OD com C20, C22 livres'","'OD with free C20, C22'"),
("'arco de dados na Órbita A (r = 1,75 km)'","'data arc in Orbit A (r = 1.75 km)'"),("'erro RMS da massa (%)'","'RMS mass error (%)'"),
("'Bennu: relacional ≈ OD até meia órbita; arco longo exige janelas'","'Bennu: relational ≈ OD up to half an orbit; long arcs need windows'"),
("'Relacional bruto'","'Relational, raw'"),("'Relacional + correção J2 fechada'","'Relational + closed-form J2 correction'"),
("'órbita da missão Psyche (altitude)'","'Psyche mission orbit (altitude)'"),
("'Psyche: o viés do achatamento cresce ao descer; a correção J2 o remove'","'Psyche: oblateness bias grows with descent; the J2 correction removes it'"),
],
'psyche_time.py':[
("'Conservador'","'Conservative'"),("'Central'","'Central'"),("'Otimista'","'Optimistic'"),("'Plano oficial'","'Official plan'"),
('d na Orbit D','d in Orbit D'),("'psyche_time.json'","'psyche_time_en.json'"),
("'dias desde o início da Orbit A (ago/2029)'","'days since start of Orbit A (Aug 2029)'"),
("'Psyche: dias de Orbit A e de margens de transferência realocados para a Orbit D (mapeamento elementar, 75 km)'","'Psyche: Orbit A days and transfer margins reallocated to Orbit D (elemental mapping, 75 km)'"),
],
'fig9.py':[
("'1/4 de órbita\\n(16 h)'","'1/4 orbit\\n(16 h)'"),("'1/2 órbita\\n(30 h)'","'1/2 orbit\\n(30 h)'"),("'1 órbita\\n(60 h)'","'1 orbit\\n(60 h)'"),
("'referência OSIRIS-REx: 4,892 ± 0,006 m³/s²'","'OSIRIS-REx reference: 4.892 ± 0.006 m³/s²'"),
("'relacional bruto'","'relational, raw'"),("'relacional + correções fechadas (offset de SRP, J2)'","'relational + closed-form corrections (SRP offset, J2)'"),
("'3ª lei de Kepler (exige 1 órbita)'","'Kepler 3rd law (needs 1 orbit)'"),("'GM de Bennu (m³/s²)'","'Bennu GM (m³/s²)'"),("'erro (%)'","'error (%)'"),
("'requisito da missão antes da Orbital A: 2 % (fora da escala, ±0,10 m³/s²)'","'mission requirement before Orbital A: 2 % (off scale, ±0.10 m³/s²)'"),
("'Efemérides reais da OSIRIS-REx (JPL Horizons, Orbital A, 3–7 jan 2019): GM de Bennu pelo estimador relacional'","'Real OSIRIS-REx ephemerides (JPL Horizons, Orbital A, 3–7 Jan 2019): Bennu GM from the relational estimator'"),
],
'fig10.py':[
("'faixa admissível (ajustável)'","'admissible band (adjustable)'"),("'GM estimado (ΔF_geom → 0)'","'estimated GM (ΔF_geom → 0)'"),
("'GM verdadeiro 1,601'","'true GM 1.601'"),("'modelo de forma 1,824 (+14 %)'","'shape model 1.824 (+14 %)'"),
("['0,4','0,6','1,0','1,6','2,0','3,0']","['0.4','0.6','1.0','1.6','2.0','3.0']"),
("'dias desde o início da Orbit A'","'days since start of Orbit A'"),
("'Massa sem informação prévia: faixa de densidade 1–8 g/cm³\\nse fecha em ±0,2 % no 1º dia'","'Mass with no prior: the 1–8 g/cm³ density band\\ncloses to ±0.2 % on day 1'"),
("'massa prévia fixa (+14 %)'","'fixed prior mass (+14 %)'"),("'espera massa do solo (7 dias)'","'wait for ground mass (7 days)'"),("'ROOP integrado (sem prévia)'","'integrated ROOP (no prior)'"),
("'alvo Orbit B1 (413 km de raio)'","'Orbit B1 target (413 km radius)'"),("'raio médio de uma órbita (km)'","'one-orbit mean radius (km)'"),
("'Descida A→B1: ROOP chega ~6 dias antes; massa errada\\nestaciona 60 km acima do alvo'","'A→B1 descent: ROOP arrives ~6 days earlier; wrong mass\\nparks 60 km above target'"),
],
'fig11.py':[
("'espera gravidade do solo (7 dias)'","'wait for ground gravity (7 days)'"),("'ROOP, J2 fixo do modelo de forma'","'ROOP, J2 fixed from shape model'"),("'ROOP integrado (massa + J2 em voo)'","'integrated ROOP (mass + J2 in flight)'"),
("'B1 → D (413 → 188 km), descida retrógrada'","'B1 → D (413 → 188 km), retrograde descent'"),
("'dias desde o início da transferência'","'days since start of transfer'"),("'raio médio de uma órbita (km)'","'one-orbit mean radius (km)'"),
("'alvo D'","'D target'"),("'alvo C'","'C target'"),
]}
for f,reps in R.items():
    s=open(f).read()
    for a,b in reps:
        if a not in s: print('MISSING',f,a[:50])
        s=s.replace(a,b)
    s=re.sub(r"savefig\('(fig\w+)\.png'",r"savefig('/home/claude/en/\1_en.png'",s)
    open('en_'+f,'w').write(s)
