# ROOP: a Invariância Relacional G′r = Gp aplicada à autonomia orbital e à estimação de massa de asteroides

Nov 21, 2026 · @Karlos

## Resumo

A Invariância Relacional G′r = Gp, com G′ = 3ω²/4πρ e ρ = 3M/4πr³ (Maia, 2026, seção 2.2), é uma reformulação exata da mecânica kepleriana. Ela está para Kepler–Newton assim como a mecânica lagrangiana está para Newton: não prevê nada diferente, mas encurta certas contas. Este artigo deriva dela o ROOP, um protocolo com duas ferramentas fechadas e embarcáveis, e as testa contra práticas de mercado e contra dados reais de voo.

**Correção orbital autônoma.** O comando relacional ΔF diz se, quanto e onde disparar. Em LEO a 350 km, com J2, J3, arrasto com atividade solar variável e ruído de GNSS, ele gasta 21,0 m/s/ano, o mínimo físico imposto pelo arrasto. Mantém a altitude média em 41 m RMS e elimina 49 a 116 ciclos anuais de planejamento em solo. Contra o solo semanal, economiza 21–27 % de Δv (1,3–1,8 kg/ano de hidrazína em 500 kg); contra o solo a cada 3 dias, \~1 %. Empata com o controle autônomo convencional usando \~25 % menos operações por amostra.

**Massa de asteroides.** O estimador M = ⟨h⟩²/(G·p\_geom) é linear e dispensa chute inicial e propagação.

- **Simulação.** Num corpo tipo Bennu, iguala a determinação de órbita (OD) em lote até \~1 %, com 300 a 2 000× menos processamento. Em Psyche, uma correção fechada de J2 o leva a 0,001 % na órbita A. Durante a espiral A→B1, entrega GM (< 0,03 %) e J2 (\~10 %) a bordo a partir do 4º dia, o que permite realocar de 2 a 24 dias para a órbita de mapeamento elementar. Em malha fechada e sem massa prévia, o ROOP usa a mesma grandeza ΔF para corrigir a massa e comandar o empuxo, e chega à Orbit B1 6 a 7 dias antes de uma estratégia que espera a massa do solo, com o mesmo Δv.
- **Dados reais de voo.** Aplicado à trajetória reconstruída da OSIRIS-REx em Bennu (JPL Horizons), dá GM = 4,897 m³/s² com uma órbita, contra a referência de 4,892 ± 0,006 m³/s² (+0,11 %, dentro de 1σ). Com janelas de 16 h, o erro fica entre −0,05 % e +0,30 %, enquanto a 3ª lei de Kepler precisa de 61 h e erra −0,36 %.

A reformulação é útil no sentido lagrangiano: não muda a física nem supera o mínimo físico, mas transforma cálculos que hoje dependem de elementos osculadores, propagação ou solo em contas fechadas, baratas e embarcáveis.

## Contribuições

Este artigo reivindica as contribuições abaixo, cada uma com a evidência que a sustenta.

1. **Leitura operacional da Invariância Relacional.** A relação G′r = Gp separa o tamanho da órbita (⟨G′r⟩ contra Gp) da forma (oscilação de K = G′/G ao longo da órbita), e cada parte se torna uma ferramenta de engenharia (seção 3).
2. **ΔF*, um comando de empuxo relacional completo.*\* Uma única grandeza entrega o critério de disparo, o tamanho do pulso em m/s e o ponto de disparo u\* que corrige a excentricidade sem Δv extra. Atinge o mínimo físico de Δv em simulação e empata com o estado da arte autônomo, com \~25 % menos operações por amostra (seção 4).
3. **Estimador relacional de massa M = ⟨h⟩²/(G·p\_geom).** Combina um ajuste puramente geométrico de p com h cinemático. É linear e fechado, vale para órbitas e sobrevoos hiperbólicos, iguala a OD em lote até \~1 % com 300 a 2 000× menos processamento e foi validado na trajetória real da OSIRIS-REx (seções 3.2, 5.2 e 5.6).
4. **Correções fechadas para perturbações reais.** O fator μ\_ef = μ\[1 + 1,5·J2·(R/r)²(1 − 3z²/r²)\] remove o viés do achatamento (0,001 % na órbita A de Psyche). O fator ⟨(r/ρ)³⟩ remove o deslocamento da órbita causado pela pressão de radiação solar, medido nos próprios dados (78 m na Orbital A de Bennu) (seções 5.4 e 5.6).
5. **Estimação conjunta de GM e J2 a bordo, inclusive durante espirais de propulsão elétrica.** μ\_ef é linear em μ e μ·J2, então estimativas em altitudes diferentes formam uma regressão linear. Na espiral A→B1 de Psyche, o GM fica abaixo de 0,03 % e o J2 em \~10 % a partir do 4º dia (seção 5.5). Integrada ao ΔF pela mesma grandeza (ΔF\_geom para a massa, dentro de uma faixa admissível ajustável, sem informação prévia), a malha leva a nave da Orbit A à B1 de Psyche 6 a 7 dias antes de esperar a massa do solo, com o mesmo Δv, e repete o ganho de \~7 dias nas transferências B1→D (com descida retrógrada, necessária para cruzar as ressonâncias de rotação) e D→C (seções 3.3, 5.7 e 5.8).
6. **Quantificação de valor sob premissas explícitas.** Δv, propelente, vida útil e ciclos de solo contra quatro estratégias de mercado (seção 4.3); dias realocáveis na missão Psyche e receita equivalente em mineração (seções 5.3 e 5.5).

**O que o artigo não reivindica:** física nova; Δv abaixo do mínimo imposto pelas perturbações; superioridade de desempenho sobre o controle autônomo convencional; precisão de massa superior à de uma OD completa; substituição do modelo de forma, do campo gravitacional de alto grau ou do rastreio radiométrico; e previsões econômicas. A relação h² = μp e o ajuste de cônicas a posições são clássicos; a contribuição está na combinação e nas correções fechadas acima, sujeita a uma revisão bibliográfica de anterioridade.

## 1. Introdução: o que torna uma reinterpretação útil

A mecânica lagrangiana não acrescentou nenhuma previsão à de Newton. Seu valor está em trocar vetores de força e reações de vínculo por um escalar e coordenadas generalizadas, o que encurta cálculos que na forma newtoniana são longos. Uma reinterpretação útil se julga, portanto, por três critérios: (i) reproduz exatamente as previsões da teoria original; (ii) leva a contas mais curtas, fechadas ou robustas para algum problema prático; (iii) esse ganho se mede contra os métodos que já se usam.

A formulação relacional reescreve as leis de Kepler em termos de um acoplamento local G′ = 3ω²/4πρ e de uma densidade orbital ρ = 3M/4πr³. Sua relação central, G′r = Gp, iguala uma grandeza medida ponto a ponto (ω, r) a uma grandeza da geometria da órbita (p). O ceticismo natural é que isso seria apenas h² = μp com outro nome. A seção 2 confirma que é — o critério (i) — e o restante do artigo testa os critérios (ii) e (iii) em dois problemas de interesse econômico:

1. **Correção orbital autônoma.** Quanto Δv, propelente e supervisão humana um controlador baseado na invariância economiza em relação ao que o mercado usa hoje: manutenção planejada em solo e controle autônomo convencional, como os demonstrados em voo no DEMETER (CNES) e no PRISMA (DLR)?
2. **Massa de asteroides em expedições de mineração.** Quanto tempo de campanha um estimador de massa baseado na invariância economiza, e quanto esse tempo vale em minério extra?

Todos os resultados vêm de simulação numérica reprodutível; código e premissas estão listados ao final. Onde a reformulação não ajuda, isso também é relatado (seção 6).

## 2. Formulação relacional

O acoplamento G′ é definido a partir de duas grandezas locais: a velocidade angular instantânea ω e a densidade orbital ρ.

A formulação desta seção é a de Maia (2026): a definição G′ ≡ ω²r³/M (seção 2.1, eq. 2), o invariante relacional G′r = Gp (seção 2.2, eq. 3) e o índice kepleriano K = G′/G = p/r = 1 + e cos θ (seção 2.3, eq. 4), em *The Orbital Coherence Function: A Keplerian Pre-Filter for Ensemble Orbit Determination* ([doi:10.5281/zenodo.20248439](https://doi.org/10.5281/zenodo.20248439)), base teórica única deste artigo. Como observa aquele trabalho, G′ não é uma nova constante física: sua utilidade é notacional, ao transformar uma lei de conservação global numa grandeza local. A densidade orbital ρ adiante é apenas uma reescrita equivalente de G′. Nada além disso é assumido; as ferramentas das seções 3 a 5 são derivadas aqui.

### 2.1 Densidade orbital

A densidade orbital é a massa do corpo central distribuída sobre a esfera de raio igual ao raio orbital instantâneo:

```latex
\rho = \frac{M}{\tfrac{4}{3}\pi r^{3}} = \frac{3M}{4\pi r^{3}}, \qquad r = R_E + h
```

ρ é uma parametrização gravitacional, não a densidade atmosférica que governa o arrasto. Ela decai como r⁻³: ≈ 4,2 × 10³ kg/m³ a 600 km de altitude e ≈ 19 kg/m³ em GEO.

| Regime | Altitude | ρ (kg/m³) |
| --- | --- | --- |
| LEO baixo | 200 km | ≈ 5,0 × 10³ |
| LEO | 600 km | ≈ 4,2 × 10³ |
| LEO alto | 2000 km | ≈ 2,4 × 10³ |
| MEO (GNSS) | 20 200 km | ≈ 7,6 × 10¹ |
| GEO | 35 786 km | ≈ 1,9 × 10¹ |

### 2.2 Acoplamento gravitacional local

O acoplamento local é definido como G′ ≡ ω²r³/M (Maia, 2026, eq. 2). Com a densidade orbital ρ, escreve-se de forma equivalente:

```latex
G' = \frac{3\,\omega^{2}}{4\pi\rho}
```

Substituindo ρ, volta-se à definição, que é a terceira lei escrita com o raio instantâneo:

```latex
G' = \frac{\omega^{2} r^{3}}{M} \quad\Longleftrightarrow\quad \omega^{2} = \frac{G' M}{r^{3}}
```

Na forma de ρ, G′ é lido como razão entre a dinâmica angular local e a massa contida na esfera orbital. Na forma kepleriana, é o valor que G teria se a órbita fosse circular com o raio e a velocidade angular daquele instante. As duas formas são algebricamente idênticas; a de ρ é a leitura relacional, a de ω²r³/M é a que se implementa.

### 2.3 Invariância Relacional

Com ω = h/r² (segunda lei de Kepler) e h² = GMp, onde p = a(1 − e²) é o semi-latus rectum:

```latex
G' r = \frac{\omega^{2} r^{4}}{M} = \frac{h^{2}}{M} = G\,a(1-e^{2}) = G\,p
```

O lado esquerdo é local (depende de r e ω no instante); o direito é global (constantes da cônica). Dimensionalmente, G′r tem unidade m⁴ kg⁻¹ s⁻², isto é, momento angular específico ao quadrado por unidade de massa central.

### 2.4 Índice kepleriano

O índice K = G′/G = p/r = 1 + e cos θ reproduz a equação polar da cônica. Sua média aritmética em anomalia verdadeira é exatamente 1, de modo que G é a média angular de G′ ao longo da órbita. O índice efetivo Kₑ = (r\_rms,θ/a)³ = (1 − e²)^¾ usa o raio quadrático médio angular r\_rms,θ = √(ab), resultado geométrico exato.

### 2.5 Leitura: reinterpretação, não física nova

A cadeia acima mostra que G′r = Gp equivale a h² = μp. Toda previsão dinâmica coincide com a de Newton–Kepler. O que a formulação oferece é uma relação explícita entre uma medida cinemática local (ω, r) e uma medida geométrica global (p). As seções seguintes testam se essa relação rende algo prático.

## 3. Duas ferramentas derivadas da invariância

A invariância separa naturalmente o **tamanho** da órbita (quanto G′r difere de Gp) da **forma** (como G′ oscila ao longo da órbita). Cada separação dá uma ferramenta.

### 3.1 ΔF: comando de empuxo relacional

A cada amostra de navegação, formam-se duas razões adimensionais:

```latex
Q_i = \frac{G'_i\, r_i}{G\, p_{ref}} = \frac{p_i}{p_{ref}}, \qquad K_i = \frac{G'_i}{G} = \frac{p}{r_i} = 1 + e\cos(u_i-\omega)
```

Sobre uma janela deslizante de uma órbita, a média de Q dá o desvio de tamanho e um ajuste linear de K contra o argumento de latitude u dá o vetor excentricidade:

```latex
\Delta F = -\frac{\bar v}{2}\big(\langle Q\rangle_T - 1\big)\ \ [\mathrm{m/s}], \qquad K_i \approx A + B\cos u_i + C\sin u_i \Rightarrow \vec e = (B/A,\, C/A)
```

A regra de acionamento: disparar quando ΔF > Δv\_min, com um pulso transversal de tamanho ΔF, centrado no argumento de latitude u\* = atan2(Δe\_y, Δe\_x) que também leva e⃗ ao alvo. Quando a correção de excentricidade excede o que o pulso de arrasto entrega, adiciona-se um par de pulsos opostos separados de meia órbita. A média de uma órbita remove a oscilação de curto período do J2 (±6 km em p a 350 km) e o ruído de GNSS.

O cálculo não passa por elementos keplerianos: usa |r × v|, |r| e u, três somas acumuladas e uma solução 3×3 por órbita.

### 3.2 Estimador relacional de massa

Lida ao contrário, a invariância dá a massa central:

```latex
M = \frac{\langle h \rangle^{2}}{G\, p_{geom}}, \qquad \frac{1}{r_i} \approx A + B\cos\theta_i + C\sin\theta_i \Rightarrow p_{geom} = 1/A
```

p\_geom vem de um ajuste puramente geométrico das posições no plano orbital, que não usa μ nem tempo; ⟨h⟩ vem da cinemática. A conta é linear, fechada e vale para órbitas elípticas e para sobrevoos hiperbólicos. Em corpos pequenos, uma aceleração de pressão de radiação solar (SRP) conhecida f é removida em primeira ordem antes do ajuste: r → r − f·Δt²/2 e v → v − f·Δt.

### 3.3 ROOP integrado: a mesma grandeza nas duas frentes

O estimador de massa e o comando de empuxo são duas leituras do mesmo desvio relacional. Como G′ = ω²r³/M exige uma massa, o protocolo trabalha com a massa estimada M̂. Comparar o lado cinemático com dois lados diferentes dá dois sinais:

```latex
\Delta F_{geom} = G'r - G\,p_{geom} = h^{2}\left(\frac{1}{\hat M} - \frac{1}{M}\right), \qquad \Delta F = -\frac{\bar v}{2}\left(\frac{\langle G'r\rangle}{G\,p_{ref}} - 1\right)
```

**ΔF\_geom** (cinemática contra geometria) é exatamente o erro de massa: zera se e somente se M̂ = M. O ROOP corrige M̂ até ΔF\_geom → 0 dentro de uma **faixa admissível ajustável**:

1. **Sem informação prévia**, a faixa inicial vem do tamanho do corpo e de densidades plausíveis (1 000 a 8 000 kg/m³). Para Psyche, GM de 0,41 a 3,0 × 10⁹ m³/s², um fator de 7.
2. Cada janela de 1/4 de órbita dá uma estimativa ⟨h⟩²/(G·p\_geom), corrigida para J2. A estimativa é limitada à faixa, de modo que uma janela espúria não tira M̂ do intervalo físico.
3. A faixa é recentrada na mediana das janelas recentes, com meia-largura igual ao maior entre 3σ robusto e 0,2 %, e nunca se alarga.
4. O empuxo só é liberado quando a meia-largura cai abaixo de 2 %. Isso torna a tolerância um parâmetro de missão, ajustável ao risco aceito.

**ΔF** (cinemática contra referência) comanda o empuxo com o M̂ corrente. Se a massa estiver errada, ΔF leva a nave a h² = G·M̂·p\_ref, ou seja, a uma órbita com p = p\_ref·M̂/M, errada na mesma proporção da massa. Por isso as duas frentes precisam rodar juntas.

## 4. Objetivo 1: correção orbital autônoma contra o mercado

O controlador com ΔF gasta o mínimo físico de Δv e elimina o planejamento em solo; a economia de Δv frente ao solo varia de \~1 % a 27 % conforme a cadência e a tolerância exigida.

### 4.1 Cenário e concorrentes

| Item | Valor |
| --- | --- |
| Órbita | 350 km, i = 51,6°, alvo: semieixo médio de referência e vetor excentricidade congelado (J2/J3, e\_f ≈ 8,7 × 10⁻⁴) |
| Forças | dois corpos + J2 + J3 + arrasto com atividade solar estocástica (fator log-normal, σ = 0,35, correlação 2 dias) |
| Satélite e propulsão | 500 kg, C\_D A/m = 0,0044 m²/kg, propulsor de 0,1 N, pulsos ≤ 0,1 m/s |
| Navegação | GNSS com 5 m e 5 cm/s por eixo (1σ) |
| Duração | 60 dias × 3 sementes; taxas anuais medidas em regime (dias 20–60) |

Cinco estratégias, todas com o mesmo propulsor, as mesmas forças e o mesmo ruído:

- **Solo semanal, ingênuo**: ciclo semanal com determinação de órbita perfeita, previsão de arrasto por persistência e 1 dia de latência de envio; pulso de arrasto em qualquer ponto, excentricidade corrigida por pares separados.
- **Solo semanal, inteligente**: idem, mas o pulso de arrasto é posicionado para corrigir também a excentricidade (boa prática).
- **Solo a cada 3 dias, inteligente**: a mesma boa prática com cadência mais alta.
- **Autônomo convencional**: a bordo, com semieixo e vetor excentricidade osculadores médios em uma órbita, na linha do DEMETER e do PRISMA.
- **Autônomo ΔF**: a bordo, com a mesma lógica de decisão do convencional; só o sensoriamento muda (Q e K relacionais).

### 4.2 Resultados

| Estratégia | Δv em regime (m/s/ano) | Excesso sobre o arrasto | Altitude média, RMS | Erro de excentricidade, RMS | Ciclos de solo/ano |
| --- | --- | --- | --- | --- | --- |
| Solo semanal, ingênuo | 28,7 ± 3,5 | +36 % | 460 m | 1,0 × 10⁻⁴ | 49 |
| Solo semanal, inteligente | 26,6 ± 6,0 | +26 % | 462 m | 1,2 × 10⁻⁴ | 49 |
| Solo a cada 3 dias, inteligente | 21,2 ± 1,8 | +0,9 % | 244 m | 4,3 × 10⁻⁵ | 116 |
| Autônomo convencional | 21,0 ± 1,7 | +0,1 % | 42 m | 5,4 × 10⁻⁶ | 0 |
| **Autônomo ΔF** | **21,0 ± 1,7** | **+0,1 %** | **41 m** | **5,5 × 10⁻⁶** | **0** |

*Média ± desvio entre 3 sementes de atividade solar e ruído. O arrasto em regime foi de 21,0 m/s/ano em todos os casos.*

&#91;image: Δv em regime, erro de altitude e ciclos de solo por estratégia\]

*Figura 1. Δv anual em regime (com o arrasto como mínimo físico), erro de altitude média e ciclos de solo por ano. Os dois autônomos coincidem; o solo só se aproxima do mínimo com cadência de 3 dias e 116 ciclos por ano.*

- **O mínimo físico é atingido.** Os dois autônomos gastam o Δv que o arrasto impõe, nem mais nem menos. Nenhuma lei de controle pode gastar menos; a pergunta é quanto as outras desperdiçam.
- **O desperdício do solo depende do requisito.** Com ciclo semanal e excentricidade exigida em 10⁻⁴, o solo gasta 26–36 % a mais, porque a previsão de arrasto erra e a excentricidade escapa da caixa entre ciclos, exigindo pares de correção. Com ciclo de 3 dias, o desperdício cai para \~1 %, ao custo de 2,4× mais ciclos de solo.
- **A precisão é o ganho robusto.** Em qualquer cadência de solo, o autônomo mantém a altitude 6–11× mais justa e a excentricidade 8–22× mais precisa.
- **ΔF empata com o convencional.** O desempenho de controle é o mesmo, e o ruído de medida também (88 m em p e em a). A diferença está no custo: o sensoriamento relacional usou 2,3 µs por amostra contra 3,0–3,5 µs do convencional na mesma implementação, \~25 % menos operações, por não precisar de energia (vis-viva) nem do vetor de Laplace.

### 4.3 Economia

Conversão de Δv em propelente para um satélite de 500 kg (m = m₀\[1 − exp(−Δv/I\_sp g₀)\]):

| Referência de comparação | Δv economizado | Hidrazína (I\_sp 220 s) | Propulsão elétrica (I\_sp 1500 s) | Vida útil, se o propelente limitar | Ciclos de solo eliminados |
| --- | --- | --- | --- | --- | --- |
| Solo semanal, ingênuo | 7,7 m/s/ano (27 %) | 1,77 kg/ano | 0,26 kg/ano | +37 % | 49/ano |
| Solo semanal, inteligente | 5,6 m/s/ano (21 %) | 1,28 kg/ano | 0,19 kg/ano | +27 % | 49/ano |
| Solo a cada 3 dias, inteligente | 0,24 m/s/ano (1 %) | 0,06 kg/ano | 0,01 kg/ano | +1 % | 116/ano |

Em escala, uma constelação de 100 satélites nessa altitude que hoje opere com ciclo semanal economizaria \~130 kg de hidrazína por ano e \~4 900 ciclos de planejamento, envio e calibração. A supervisão passa a ser monitoramento: no DEMETER, o acompanhamento de rotina do controle autônomo tomava alguns segundos de gráficos de telemetria por semana. O custo-hora desses ciclos varia por operador e não foi estimado aqui; o número de ciclos eliminados é a grandeza a multiplicar pelo custo de cada um.

Atribuição honesta: a economia frente ao solo é da **autonomia em malha fechada**, que um controlador convencional também entrega. O que a formulação relacional acrescenta é o caminho mais curto até ela: uma única grandeza (G′) dá tamanho, forma e ponto de disparo, com menos operações.

## 5. Objetivo 2: massa de asteroides e retorno de expedições

O estimador relacional atinge a mesma precisão de uma OD completa com \~300× menos processamento e sem chute inicial; o que ele economiza é o tempo de ida e volta ao solo, e esse tempo se converte em dias de mineração.

### 5.1 Como o mercado faz hoje

A massa do corpo-alvo condiciona a inserção em órbita e as descidas. Na OSIRIS-REx, a fase Preliminary Survey (cinco sobrevoos hiperbólicos a \~7,25 km de Bennu) tinha como requisito uma massa com erro < 1–2 % antes da inserção na órbita Orbital A, em 31/12/2018. A estimativa prévia em solo tinha \~11,5 % de incerteza; a primeira metade do primeiro sobrevoo já a reduziu a \~0,7 %, e o fim da fase a \~0,14 %. Cada solução de OD saiu cerca de 24 h depois do início do trecho correspondente. Na Hayabusa2, a massa de Ryugu veio de uma operação dedicada de descida em queda livre, em 6–7/08/2018, cerca de seis semanas depois da chegada à posição de estacionamento.

Dois fatos orientam o teste: a precisão é limitada pelos dados, não pelo algoritmo; e o tempo é consumido pelo ciclo dados → solo → decisão → envio.

### 5.2 Teste em corpo tipo Bennu

GM = 4,89 m³/s²; SRP de 1 × 10⁻⁷ m/s² sobre a nave (comparável à gravidade a 7 km); ruído de 2 m em posição e 0,1 mm/s em velocidade; 60 medidas por arco; 40 Monte Carlo por ponto. Concorrente: OD em lote com μ e seis estados livres, propagação numérica e SRP modelada.

| Geometria | Arco | Relacional, SRP compensada (1σ) | OD em lote (1σ) | Relacional sem compensar SRP |
| --- | --- | --- | --- | --- |
| Órbita a 1,75 km (período 58 h) | 4 h | 2,0 % | 2,0 % | 5,3 % |
| Órbita a 1,75 km | 8 h | 0,55 % | 0,57 % | 4,7 % |
| Órbita a 1,75 km | 16 h | 0,35 % | 0,19 % | 4,6 % |
| Sobrevoo a 7,25 km, v∞ = 0,1 m/s | 12 h | 3,7 % | 3,7 % | falha |
| Sobrevoo a 7,25 km | 24 h | 0,99 % | 1,00 % | falha |
| Sobrevoo a 7,25 km | 48 h | 0,31 % (viés −0,46 %) | 0,32 % | falha |

&#91;image: Erro da massa contra duração do arco, órbita e sobrevoo\]

*Figura 2. Erro da massa (1σ) contra a duração do arco de dados. As curvas do estimador relacional com SRP compensada e da OD em lote se sobrepõem até \~1 %; sem compensação, o relacional estaciona em \~5 % em órbita.*

- **Mesma precisão da OD** até \~1 %, nos dois tipos de trajetória. Abaixo disso, a compensação de SRP de primeira ordem começa a limitar (viés de −0,46 % em 48 h de sobrevoo; 0,35 % contra 0,19 % em 16 h de órbita).
- **A SRP não pode ser ignorada** em corpos de centenas de metros: sem compensação, o estimador erra 5 % em órbita e falha no sobrevoo.
- **Custo.** O estimador relacional levou 0,29 ms por estimativa, contra \~80 ms da OD em lote na mesma máquina (\~300×), sem precisar de modelo de propagação. Neste teste simples, a OD também convergiu com chute inicial errado em até 200 %, então a vantagem demonstrada é de custo e simplicidade, não de robustez.
- **Corpos grandes e campo irregular.** Num corpo do porte de Psyche (SRP desprezível), o erro cai abaixo de 0,01 % em 1/4 de órbita. Com J2 do corpo, o estimador e a OD de massa pontual sofrem o mesmo viés (até +2 % para J2 = 0,1 a 300 km); removê-lo exige uma OD com J2 livre.

### 5.3 Tempo economizado e valor da expedição

Como a precisão é a mesma, o arco de dados necessário é o mesmo. O ganho é tirar o solo do caminho crítico: a nave estima a massa a bordo assim que o arco termina e pode seguir para a próxima fase. Modelo:

```latex
D = N_{dec}\cdot L, \qquad \Delta R = D\cdot P\cdot V
```

com N\_dec decisões dependentes da massa (inserção, rebaixamentos de órbita, aproximação para ancoragem; 3 a 5), L a latência de cada ciclo de solo (1 dia na OSIRIS-REx, até 3 dias com acesso compartilhado a antenas), P a produção diária e V o valor por kg. Isso dá **D = 3 a 15 dias**. Com a fase de mineração de 100 dias adotada por Hein et al., são 3 a 15 % a mais de operação, desde que a janela de partida seja fixa e os dias poupados virem dias de mineração.

Produtividade e preços: a vazão otimista de Hein et al. para uma nave de 150 kg dá \~149 kg/dia de água (5 % de rendimento) ou \~45 kg/dia de platina (razão 10⁻⁵); o cenário conservador usa 10 % disso. Água vale US$ 20 mil/kg entregue em órbita cislunar na premissa de Hein et al.; platina, \~US$ 43,7 mil/kg no preço à vista de 19/11/2026.

| Produto | Cenário | Produção (kg/dia) | Valor por dia | D = 3 dias | D = 9 dias | D = 15 dias |
| --- | --- | --- | --- | --- | --- | --- |
| Água | otimista | 149 | US$ 2,98 mi | US$ 8,9 mi | US$ 26,8 mi | US$ 44,7 mi |
| Água | conservador | 15 | US$ 0,30 mi | US$ 0,9 mi | US$ 2,7 mi | US$ 4,5 mi |
| Platina | otimista | 45 | US$ 1,98 mi | US$ 5,9 mi | US$ 17,8 mi | US$ 29,7 mi |
| Platina | conservador | 4,5 | US$ 0,20 mi | US$ 0,6 mi | US$ 1,8 mi | US$ 3,0 mi |

**Propelente.** A economia de propelente neste objetivo é desprezível: as manobras de proximidade em corpos pequenos são de cm/s a m/s, pequenas diante do Δv de transferência. O valor econômico do estimador está no tempo, não no combustível. A receita também não é lucro: não foram descontados custos de operação nem de retorno do material.

### 5.4 Casos reais: Bennu e Psyche com campo irregular e rotação

Com parâmetros publicados e o campo gravitacional de grau 2 girando com o corpo, o estimador relacional atinge < 1 % em Bennu em poucas horas. Em Psyche, uma correção fechada de J2 o leva a 0,003 %, a mesma qualidade de uma OD completa com C20 e C22 livres, com 300 a 2 000× menos processamento.

| Parâmetro | Bennu | Psyche |
| --- | --- | --- |
| GM | 4,892 m³/s² (OSIRIS-REx) | 1,601 km³/s² (astrometria, 2024) |
| Raio de referência | 244 m | 113,2 km |
| C20, C22 (não normalizados) | −0,0334; 0,0028 (elipsoide 252,8 × 246,2 × 228,7 m, densidade uniforme) | −0,1027; 0,0190 (momentos de inércia do modelo de forma) |
| Período de rotação | 4,30 h | 4,196 h |
| Trajetórias | sobrevoo a 7,25 km; órbitas a 1,75 km e 1,0 km | órbitas A (700 km), B (300 km), C (190 km), D (75 km) |
| SRP na nave (premissa) | 1 × 10⁻⁷ m/s² | 2,6 × 10⁻⁸ m/s² |
| Ruído de navegação (premissa) | 2 m; 0,1 mm/s | 50 m; 0,1 mm/s |
| Chute inicial da OD | 5,2 m³/s² (pré-chegada) | 1,824 km³/s² (modelo de forma) |

*A fonte dos coeficientes de Psyche imprime C22 = 0,1904; os momentos de inércia do próprio artigo dão 0,0190, valor usado aqui. 25 Monte Carlo por ponto, inclinações de 70° a 110° em relação ao eixo de rotação.*

**Bennu (erro médio ± 1σ):**

| Trajetória | Arco | Relacional | OD massa pontual | OD com C20, C22 livres |
| --- | --- | --- | --- | --- |
| Sobrevoo a 7,25 km | 24 h | +0,04 ± 0,96 % | +0,03 ± 0,89 % | +0,71 ± 4,06 % |
| Sobrevoo a 7,25 km | 48 h | −0,03 ± 0,51 % | +0,01 ± 0,27 % | +0,17 ± 0,69 % |
| Órbita a 1,75 km | 7,2 h (1/8) | +0,09 ± 0,60 % | +0,13 ± 0,64 % | +0,88 ± 2,15 % |
| Órbita a 1,75 km | 14,4 h (1/4) | −0,01 ± 0,29 % | +0,02 ± 0,20 % | −0,04 ± 0,64 % |
| Órbita a 1,75 km | 57,8 h (1) | +0,05 ± 3,51 % | −0,06 ± 0,05 % | 0,00 ± 0,06 % |
| Órbita a 1,0 km | 6,2 h (1/4) | −0,07 ± 0,28 % | −0,02 ± 0,36 % | −0,15 ± 0,94 % |

- Em Bennu o campo irregular é pequeno (C20 = −0,033 a 1,75 km dá viés < 0,15 %), e o estimador relacional acompanha a OD de massa pontual até meia órbita: < 1 % com 7 h em órbita ou 24 h de sobrevoo.
- Em arcos longos ele degrada (3,5 % numa órbita completa), porque a compensação de SRP de primeira ordem deixa de valer quando ½ f Δt² chega a centenas de metros. A regra prática é estimar em janelas de 1/4 a 1/2 órbita.
- Com arcos curtos, a OD com C20 e C22 livres fica pior que a de massa pontual: parâmetros demais para os dados.

**Psyche (órbita completa, erro médio ± 1σ):**

| Órbita | Período | Relacional bruto | Relacional + correção J2 | OD massa pontual | OD com C20, C22 livres |
| --- | --- | --- | --- | --- | --- |
| A, 700 km | 32,0 h | −0,13 ± 0,02 % | +0,001 ± 0,003 % | −0,14 ± 0,02 % | 0,000 ± 0,003 % |
| B, 300 km | 11,6 h | −0,49 ± 0,07 % | +0,007 ± 0,012 % | −0,54 ± 0,07 % | −0,001 ± 0,005 % |
| C, 190 km | 7,3 h | −0,87 ± 0,26 % | +0,04 ± 0,20 % | −1,18 ± 0,43 % | −0,001 ± 0,007 % |
| D, 75 km | 3,6 h | −2,22 ± 1,02 % | +0,07 ± 0,82 % | −2,88 ± 1,22 % | −0,001 ± 0,013 % |

- **O achatamento domina em Psyche.** Com C20 = −0,103, o estimador relacional e a OD de massa pontual medem um μ efetivo, e o viés cresce de −0,13 % na órbita A para −2,2 % na D.
- **A correção é fechada.** Para órbita quase circular com inclinação i em relação ao eixo de rotação:

```latex
\mu_{ef} = \mu\left[1 + \tfrac{3}{2}J_2\left(\tfrac{R}{r}\right)^{2}\left(1 - \tfrac{3}{2}\sin^{2} i\right)\right], \qquad J_2 = -C_{20}
```

Dividir a estimativa relacional por esse fator, com J2 do modelo de forma, leva o erro a 0,001 % na órbita A e 0,007 % na B. Com J2 errado em 20 %, o erro fica em 0,03 % e 0,11 %. Nas órbitas C e D sobra uma dispersão de 0,2–0,8 % vinda do C22 em rotação, que a correção não cobre.

- **Duas altitudes resolvem massa e J2 sem OD.** Como μ\_ef é linear em μ e μ·J2, as estimativas relacionais em duas órbitas formam um sistema 2×2. Com A e C: massa com −0,005 ± 0,03 % e J2 com −6 ± 27 %.
- **Custo.** O estimador relacional levou 0,3 ms por estimativa; a OD com C20 e C22 livres, 100 a 600 ms por estimativa na mesma máquina.
- **Contexto.** A massa de Psyche já é conhecida da Terra a \~1 % por astrometria (GM = 1,601 ± 0,017 km³/s²), e o modelo de forma com densidade uniforme sugere \~14 % a mais. A missão pretende chegar a GM melhor que 0,001 km³/s² com rastreio radiométrico. O estimador relacional com correção J2 já entrega essa ordem de precisão (0,003 %, ou \~5 × 10⁻⁵ km³/s²) a bordo, na primeira órbita.

&#91;image: Erro da massa em Bennu por arco e em Psyche por órbita\]

*Figura 3. Esquerda: Bennu, órbita a 1,75 km; o relacional acompanha a OD até meia órbita. Direita: Psyche, órbita completa em cada fase da missão; o viés do achatamento cresce ao descer e a correção J2 fechada o remove nas órbitas A e B.*

### 5.5 Psyche: do cálculo de M ao tempo de exploração mineral

Na Psyche, a exploração mineral propriamente dita é a Orbit D: o mapeamento elementar a 75 km. Estimar GM e J2 a bordo, continuamente, inclusive durante as espirais de propulsão elétrica, pode realocar de 2 a 24 dias para essa órbita, ou seja, +1 a +15 % sobre os 163 dias planejados. O cenário central dá +12 dias.

**Onde a massa é portão.** No plano oficial, a Orbit A dura 56 dias, contra um mínimo de 40, e serve para caracterizar gravidade e forma. Os dados de gravidade e forma da Orbit A são necessários para projetar o empuxo até a Orbit B, e cada mudança de altitude é limitada pelo conhecimento de gravidade disponível. Cada início de campanha recebe horas ou dias de margem para a incerteza do tempo de transferência. As transferências somam 228 dias (17 + 98 + 90 + 23) dos \~740 dias de fase orbital.

**Teste: estimador relacional durante a espiral Orbit A → B1.** Espiral de 811 km para 413 km de raio em 17 dias, como no plano; campo de grau 2 em rotação, SRP e ruído de 50 m e 0,1 mm/s. Em cada janela de 1/4 de órbita, as acelerações conhecidas (empuxo comandado e SRP) são removidas por dupla integração e o estimador dá μ\_ef(r). Uma regressão linear de μ\_ef contra g = 1,5⟨(R/r)²(1 − 3z²/r²)⟩ dá GM e J2. 8 realizações:

| Dia da espiral | Altitude | Erro de GM | Erro de J2 |
| --- | --- | --- | --- |
| 1 | 661 km | −0,009 ± 0,026 % | −6 ± 16 % |
| 2 | 638 km | +0,001 ± 0,013 % | +1 ± 8 % |
| 4 | 569 km | +0,013 ± 0,006 % | +7 ± 3 % |
| 8 | 459 km | +0,018 ± 0,002 % | +8 ± 1 % |
| 17 | 303 km | +0,026 ± 0,007 % | +10 ± 4 % |

O GM fica abaixo de 0,03 % desde o primeiro dia, e o J2 se estabiliza em \~10 % a partir do 4º dia (o viés restante vem do C22 em rotação). Erro de 1 % ou 3 % no conhecimento do empuxo não alterou os resultados: a aceleração de empuxo (\~1 × 10⁻⁵ m/s²) é duas ordens menor que a gravidade nessas altitudes. Na prática, a nave carrega um modelo de grau 2 atualizado a bordo ao longo de cada transferência, sem esperar o modelo do solo.

**Tradução em dias.** Os dias poupados vêm de duas fontes: a folga da Orbit A sobre o mínimo (até 16 dias, se a gravidade for o item limitante) e as margens de início das quatro campanhas pós-transferência (0,5 a 2 dias cada, premissa). O total é realocado para a Orbit D:

| Cenário | Corte na Orbit A | Margens de transferência | Dias a mais na Orbit D | Sobre os 163 dias planejados | Sobre o mínimo de 100 dias | Órbitas extras a 75 km | Incerteza de abundância elementar |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Conservador | 0 d | 4 × 0,5 d | +2 d | +1,2 % | +2 % | +14 | −0,6 % |
| Central | 8 d | 4 × 1 d | +12 d | +7,4 % | +12 % | +82 | −3,5 % |
| Otimista | 16 d | 4 × 2 d | +24 d | +14,7 % | +24 % | +164 | −6,6 % |

*Órbitas extras pela taxa da Orbit D (684 órbitas em 100 dias). A incerteza de contagem do espectrômetro cai com 1/√T.*

&#91;image: Linha do tempo da Psyche: plano oficial e três cenários com dias realocados para a Orbit D\]

*Figura 4. Fases da Psyche no plano oficial e nos três cenários; em verde escuro, os dias realocados para a Orbit D.*

**Limites honestos.**

- **O que o estimador não substitui.** Ele não substitui o modelo de forma (imagens) nem o campo gravitacional de grau 10 da Orbit C, que exige pelo menos 35 dias de rastreio radiométrico; esses tempos não encolhem.
- **O corte na Orbit A.** Ele só existe se a gravidade, e não a cobertura de imagens, for o que segura a fase. O cenário otimista supõe que é.
- **A margem por campanha.** A fonte não dá o valor exato, só "horas ou dias", por isso a faixa de 0,5 a 2 dias.
- **O fim da missão.** A data de término é fixa; os dias ganhos são realocação, não extensão.
- **Em mineração comercial.** Com as produtividades da seção 5.3 (US$ 0,2 a 3 milhões por dia), 2 a 24 dias equivalem a US$ 0,4 a 72 milhões por expedição, com a mesma ressalva de que isso é receita bruta.

### 5.6 Validação com dados reais: a trajetória da OSIRIS-REx em Bennu

Aplicado à trajetória real da OSIRIS-REx na órbita Orbital A, o estimador relacional recupera o GM de Bennu com +0,11 % de erro a partir de uma órbita (dentro de 1σ da referência), e com −0,05 % a +0,30 % a partir de 16 h de dados. Nos mesmos dados, a 3ª lei de Kepler precisa de 61 h e erra −0,36 %.

**Dados.** 49 vetores de estado (posição e velocidade) da nave em relação ao centro de Bennu, obtidos do sistema JPL Horizons (alvo −64, centro @2101955, eclíptica J2000, TDB). A cadência é de 2 h, de 3/01/2019 00:00 a 7/01/2019 12:00, com uma lacuna de 14 h; isso dá um segmento contínuo de 72 h e outro de 22 h, \~1,8 órbita. Checagens de integridade: as velocidades batem com as diferenças finitas das posições, e |h| se conserva a 0,16 %. A órbita reconstruída tem período de 61,1 h, semieixo de 1 814 m e e = 0,127. Referência: GM = 4,892 ± 0,006 m³/s².

**Método, todo fechado.**

1. Plano orbital por mínimos quadrados, sem forçá-lo a passar pelo centro de Bennu.
2. Estimativa relacional no plano: M = ⟨h⟩²/(G·p\_geom), com p\_geom do ajuste de 1/ρ.
3. Correção do deslocamento fora do plano: a gravidade radial no plano vale μρ/r³, então μ = μ\_ajuste·⟨(r/ρ)³⟩.
4. Correção de J2 com C20 = −0,0334 (elipsoide ajustado de Bennu) e eixo de rotação próximo do polo da eclíptica.

O próprio ajuste mostrou que a Orbital A é uma órbita de terminador: a normal do plano fica a 0,8° da direção do Sol. A órbita está deslocada 78 m do centro de Bennu, o que implica uma aceleração de pressão de radiação de \~6 × 10⁻⁸ m/s². Esse é um parâmetro físico da nave, medido aqui sem modelo de SRP.

| Arco de dados | Janelas | Relacional bruto | Relacional + correções fechadas |
| --- | --- | --- | --- |
| 1/4 de órbita (16 h) | 9 | 4,890 ± 0,013 m³/s² (−0,05 % ± 0,27 %) | 4,906 ± 0,011 m³/s² (+0,30 % ± 0,22 %) |
| 1/2 órbita (30 h) | 3 | 4,889 ± 0,005 m³/s² (−0,07 % ± 0,09 %) | 4,903 ± 0,006 m³/s² (+0,23 % ± 0,12 %) |
| 1 órbita (60 h) | 1 | 4,882 m³/s² (−0,20 %) | **4,897 m³/s² (+0,11 %)** |
| Todos os 49 pontos (108 h, com a lacuna) | 1 | 4,879 m³/s² (−0,28 %) | 4,898 m³/s² (+0,12 %) |
| 3ª lei de Kepler, 4π²a³/T² | 1 órbita (61 h) | 4,875 m³/s² (−0,36 %) | — |

&#91;image: GM de Bennu estimado em cada janela dos dados reais da OSIRIS-REx contra a referência\]

*Figura 5. GM de Bennu estimado em cada janela da trajetória real da OSIRIS-REx (Orbital A, 3–7/01/2019). Faixa cinza: referência 4,892 ± 0,006 m³/s². O requisito da missão antes da inserção na Orbital A era 2 % (±0,10 m³/s²), fora da escala.*

**Leitura.**

- Em todas as janelas o erro fica entre −0,5 % e +0,6 %, de 3 a 40 vezes abaixo do requisito de 2 % que a missão tinha antes de entrar na Orbital A.
- Com uma órbita completa, as duas correções fechadas levam o estimador para dentro da incerteza da referência.
- Em arcos curtos, a correção de deslocamento passa do ponto (+0,2 a +0,3 %), porque o ajuste de cônica num arco parcial absorve parte dos termos não modelados (SRP no plano, C22 em rotação). O bruto fica mais perto da referência nesses arcos, mas por compensação de erros.
- Contra a solução rápida clássica, a 3ª lei de Kepler, o estimador relacional precisa de 1/4 do arco para a mesma ordem de erro e não depende de medir o período.

**Alcance do teste.** As efemérides do Horizons são a trajetória reconstruída pela navegação da missão, que usou um valor de GM e modelos de força. Por isso, este teste verifica que o estimador recupera, por geometria e cinemática apenas, a dinâmica contida numa trajetória real voada, com SRP, campo irregular e excentricidade reais. Não é uma medição independente do GM. O passo seguinte para independência total é aplicar o estimador aos dados brutos de navegação óptica e radiométrica arquivados no PDS.

### 5.7 ROOP integrado na Psyche: descida sem massa prévia

Com as duas frentes em malha fechada, a nave desce da Orbit A para a B1 sem nenhuma massa prévia e chega 6 a 7 dias antes de uma estratégia que espera a massa do solo, com o mesmo Δv e com raio final igual ou mais preciso. Usar a massa do modelo de forma (+14 %) sem correção leva a uma órbita 59–63 km acima do alvo.

**Cenário.** Psyche com campo de grau 2 em rotação, SRP e ruído de navegação de 50 m e 0,1 mm/s. O empuxo contínuo de 1,22 × 10⁻⁵ m/s² reproduz a transferência nominal de 17 dias. A descida vai de 811 km para 413 km de raio, seguida de manutenção com pulsos posicionados em u\*, num horizonte de 30 dias com 3 sementes. A chegada é contada pela trajetória verdadeira: média de uma órbita a até 5 km do alvo. Três estratégias:

- **Massa prévia fixa:** GM do modelo de forma (1,824 × 10⁹ m³/s²), descida imediata.
- **Espera a massa do solo:** fica na Orbit A até receber uma solução de solo com erro de 0,1 %, após 7 dias (premissa de latência), e então desce.
- **ROOP integrado:** sem massa prévia; ΔF\_geom corrige a massa dentro da faixa ajustável e ΔF comanda a descida assim que a faixa cai abaixo de ±2 %.

| Estratégia | Chegada à B1 | Δv até a chegada | Erro final de raio | Erro final de GM |
| --- | --- | --- | --- | --- |
| Massa prévia fixa (+14 %) | não chega (estaciona acima do alvo) | — | +59 a +63 km | +13,9 % |
| Espera a massa do solo (7 dias) | dia 24,2–25,2 | 17,8–18,9 m/s | +3,3 a +3,5 km | 0,1 % |
| **ROOP integrado, sem prévia** | **dia 18,1–18,2** | **17,9–18,0 m/s** | **+2,5 a +3,0 km** | **+0,10 a +0,17 %** |

&#91;image: Faixa admissível de massa se fechando no primeiro dia e descida A→B1 das três estratégias\]

*Figura 6. Esquerda: GM estimado e faixa admissível do ROOP sem massa prévia; a faixa de densidade 1–8 g/cm³ se fecha em ±0,2 % no primeiro dia. Direita: raio médio de uma órbita na descida A→B1 (semente 0).*

**Leitura.**

- A primeira janela de \~8 h na Orbit A já dá GM com 0,03 % de erro. Ao fim do primeiro dia a faixa está em ±0,2 % e a descida é liberada.
- O ganho de 6 a 7 dias é a latência do solo menos o dia de convergência do ROOP. Com outra latência, o ganho muda na mesma medida.
- Durante a descida, a mesma malha continua atualizando GM e J2 (seção 5.5). Por isso o raio final fica tão ou mais preciso que o da estratégia que recebeu a massa do solo, sem gastar mais Δv para chegar.
- **Aplicação à missão.** Na Psyche, os portões que dependem de gravidade são a saída da Orbit A e a descida a partir da B1. Com o ROOP integrado em cada um, a espera pelo solo deixa de estar no caminho crítico, consistente com o cenário central de +12 dias na órbita de mapeamento elementar (seção 5.5). A Orbit A continua tendo de cumprir seu mínimo por imageamento e forma; o que sai do caminho crítico é a espera pela massa.
- O controle de manutenção na B1 é o do banco de testes, não otimizado (0,6–1,0 m/s por dia contra o C22 a 300 km). Os resultados de chegada não dependem dele.

### 5.8 A malha nas transferências seguintes: B1 → D e D → C

Repetida nas duas transferências seguintes da Psyche, a malha integrada chega à órbita-alvo 7 dias antes, em média, de uma estratégia que espera a gravidade do solo, com Δv equivalente e erro de raio abaixo de 1 km. Abaixo da B1, isso só é possível com uma descida retrógrada em relação à rotação de Psyche.

**O que muda abaixo da B1.**

- O raio síncrono de Psyche, onde o período orbital iguala a rotação de 4,2 h, fica em \~210 km, e a ressonância 3:2 em \~275 km. A Orbit D (188 km) está abaixo das duas.
- No campo de grau 2 em rotação, uma descida polar atravessa essas ressonâncias e tem a excentricidade bombeada. Sem amortecimento da excentricidade, a nave chega à D no dia 3, mas nas semanas seguintes o raio mínimo (amostrado a cada hora) cai a 60, 99 e 123 km em três sementes; o raio médio de Psyche é 113 km, então duas das três terminam em impacto. Com pulsos de amortecimento em u\*, a órbita fica segura, mas estaciona em \~240 km e não chega à D em 25 dias.
- A estimativa de massa, porém, fica abaixo de 0,8 % de erro em todas as descidas polares. A falha é de guiagem, não de estimação.
- Uma descida retrógrada em relação à rotação resolve o problema, coerente com a estabilidade conhecida de órbitas retrógradas em torno de corpos alongados. A 120° e a 150° de inclinação, a nave chega à D com raio médio de uma órbita entre 186 e 191 km; a 150°, o raio mínimo (amostrado a cada hora) fica em 178 km ou mais, com ou sem amortecimento. Perto do equador (170°), o J2 deixa de ser separável do C22 e a estimativa de massa degrada para −6,5 %, um diagnóstico que a própria malha fornece.
- **Alcance.** Essa instabilidade é deste modelo e deste controlador (só campo de grau 2, descida de 3 dias com empuxo alto, sem controle de excentricidade durante a descida). Não é uma previsão sobre a missão real, que usa órbitas quase polares na D, transferências de 90 a 98 dias e o campo gravitacional completo (`psyche_polar_check.py`).

**Cenário.** Inclinação de 150° (premissa de projeto validada acima). Autoridade de empuxo de 1 × 10⁻⁴ m/s², da classe de um propulsor SPT-140 (\~0,27 N em \~2 700 kg); as durações reais das transferências (98 e 90 dias) incluem mudança de plano e ciência, não modeladas aqui. Massa herdada da fase anterior (+0,15 %); J2 inicial do modelo de forma com 20 % de erro. Excentricidade medida sobre a última órbita completa, com amortecimento acima de 0,04. A chegada é contada com raio médio a até 1 % do alvo. São três estratégias com 3 sementes cada:

- **Espera a gravidade do solo:** 7 dias na órbita de partida até receber GM com 0,1 % e J2 com 1 %.
- **ROOP com J2 fixo:** parte imediatamente; massa em voo, J2 do modelo de forma (ablação).
- **ROOP integrado:** parte imediatamente; massa e J2 em voo, por ΔF\_geom e regressão entre altitudes.

| Perna | Estratégia | Chegada (dias) | Δv até a chegada | Erro final de raio | Erro de GM | Erro de J2 |
| --- | --- | --- | --- | --- | --- | --- |
| B1 → D | Espera do solo | 14,5–17,4 | 33–37 m/s | +0,2 a +0,5 km | 0,1 % | 1 % |
| B1 → D | ROOP, J2 fixo | 6,2–12,6 | 32–38 m/s | −0,2 a +0,4 km | −0,5 a −1,1 % | 20 % |
| B1 → D | **ROOP integrado** | **6,2–13,7 (média 8,7)** | **32–39 m/s** | **−0,3 a −0,2 km** | **−0,3 a −0,4 %** | **+10 a +12 %** |
| D → C | Espera do solo | 11,0–11,2 | 22–23 m/s | −0,2 a +1,2 km | 0,1 % | 1 % |
| D → C | ROOP, J2 fixo | 4,0–4,2 | 21–23 m/s | −0,5 a +0,4 km | −0,2 a −0,3 % | 20 % |
| D → C | **ROOP integrado** | **4,0–4,2** | **21–22 m/s** | **−0,2 a +0,9 km** | **−0,1 a −0,2 %** | **+10 a +11 %** |

&#91;image: Raio médio nas transferências B1→D e D→C para as três estratégias\]

*Figura 7. Raio médio de uma órbita nas transferências B1→D (descida retrógrada) e D→C, semente 0. As duas variantes do ROOP partem imediatamente e se sobrepõem; a estratégia de solo parte após 7 dias.*

**Leitura.**

- **D → C.** O ganho é de 7,0 dias em todas as sementes, exatamente a latência do solo, com o mesmo Δv.
- **B1 → D.** O ganho médio é de 7,3 dias, mas varia de 0,9 a 11,2. Quando a passagem pelas ressonâncias exige episódios de amortecimento, a descida se alonga, e isso independe de quem estima a gravidade.
- **J2 em voo.** Estimar o J2 durante a descida reduz o erro de massa nas órbitas baixas de −0,5 a −1,1 % para −0,3 a −0,4 %. O J2 estimado carrega +10 % de viés, por aliasing do C22 em rotação.
- **Missão completa.** Somando as três pernas (A→B1, B1→D e D→C), o ROOP integrado tira cerca de 20 dias de espera pelo solo do caminho crítico, com a premissa de 7 dias de latência em cada portão. Isso fica dentro da faixa de 2 a 24 dias realocáveis estimada na seção 5.5, perto do cenário otimista.

## 6. Onde a reformulação não ajuda

Os testes negativos delimitam o uso correto da invariância; todos foram reproduzidos nas mesmas simulações.

| Uso testado | Resultado | Lição |
| --- | --- | --- |
| ΔF = G′r − Gp com p tirado de elementos osculadores do mesmo estado | ΔF ≡ 0 para qualquer estado ( | ΔF/Gp |
| ΔF avaliado a cada instante contra p de referência | Gasta \~4× o Δv do arrasto e eleva a órbita média em 6 km | O J2 faz p oscilar ±6 km a cada meia órbita; só a média de uma órbita é sinal de erro |
| Limiar fixo em unidades de G′r, abaixo do ruído | Dispara por ruído de GNSS (88 m em p) | O limiar deve ser expresso em m/s de Δv e ficar acima do ruído filtrado |
| Estimador de massa por média temporal de ω²r³/G | Viés de −e² com amostragem uniforme no tempo | A média temporal de K = 1 + e cos θ vale 1 − e²; o ajuste geométrico de p elimina o viés |
| Estimador relacional sem compensar SRP em corpo pequeno | Erro de \~5 % em órbita; falha em sobrevoo | Acelerações não gravitacionais conhecidas devem ser removidas antes do ajuste |
| Ganho de Δv sobre o mínimo físico | Nenhum controlador ficou abaixo do arrasto | A economia de Δv existe apenas contra estratégias que desperdiçam |

## 7. Conclusões

### As leis de Kepler relidas pela Invariância Relacional

As leis de Kepler são enunciados globais: a órbita inteira é uma elipse, áreas iguais são varridas em tempos iguais, e T² ∝ a³ relaciona o período ao tamanho de toda a órbita. A Invariância Relacional G′r = Gp (Maia, 2026, seção 2.2) é *relacional* no sentido literal: relaciona o global ao local e instantâneo. O lado esquerdo é medido em um único ponto (r, ω); o direito é uma constante da órbita inteira (p). Como a igualdade vale em todo ponto, o conteúdo das leis globais fica disponível a cada instante. Não é preciso derivar a trajetória, integrar equações de movimento nem esperar que a órbita se complete.

| Lei de Kepler | Forma global | Forma relacional, local e instantânea | Uso no ROOP |
| --- | --- | --- | --- |
| 1ª (órbitas cônicas) | a órbita inteira é uma elipse de parâmetros p, e | K = G′/G = p/r = 1 + e cos θ em cada ponto | forma da órbita e ponto de disparo u\* |
| 2ª (áreas iguais) | h constante ao longo da órbita | G′r = Gp em cada ponto | tamanho da órbita (ΔF) e massa (ΔF\_geom) |
| 3ª (T² ∝ a³) | GM = 4π²a³/T², uma órbita inteira | G′ ≡ ω²r³/M com r e ω instantâneos | acoplamento local que liga as duas linhas acima |

**A analogia lagrangiana.** A mecânica lagrangiana é a mecânica de Newton sem a necessidade de escrever explicitamente F = ma: a dinâmica sai de um escalar, e problemas com vínculos ficam curtos. Da mesma forma, a Invariância Relacional é a mecânica de Kepler sem a necessidade de explicitar as leis em sua forma global: uma única grandeza local, G′, carrega as três leis ponto a ponto. Em nenhum dos dois casos há física nova; nos dois, o ganho está no que se torna curto de calcular.

**Por que a 3ª lei sozinha não basta.** Ela entrega um número (GM) por órbita completa e nada mais:

- Precisa do período e do semieixo maior; nos dados reais da OSIRIS-REx, levou 61 h e errou −0,36 %. O estimador relacional chegou à mesma ordem de erro com 16 h.
- Não informa onde disparar, não separa tamanho de forma e não se aplica a sobrevoos hiperbólicos, onde não há período.

Os resultados deste artigo vêm de reunir as três leis numa única grandeza local. A massa sai de h² contra p, ou seja, da 2ª lei contra a 1ª. O tamanho do pulso sai de ⟨G′r⟩ (2ª lei), e o ponto de disparo sai da oscilação de K (1ª lei). A 3ª lei, escrita de forma instantânea em G′, é o elo que permite trabalhar com frações de órbita, tanto em órbita fechada quanto em sobrevoo.

### Síntese dos resultados

A Invariância Relacional da seção 2.2 de Maia (2026) passa nos três critérios de uma reinterpretação útil. Reproduz exatamente a dinâmica kepleriana. Gera duas ferramentas fechadas e embarcáveis. E, medida contra o mercado e contra dados reais de voo, entrega ganhos quantificáveis, todos vindos de levar a decisão para bordo com contas curtas e nenhum de física nova.

1. **Correção orbital.**
   - ΔF atinge o mínimo físico de Δv (21,0 m/s/ano a 350 km), com altitude 6–11× mais justa e zero ciclos de solo.
   - Contra o solo semanal, economiza 21–27 % de Δv (1,3–1,8 kg/ano de hidrazína em 500 kg; +27–37 % de vida útil quando o propelente limita). Contra o solo a cada 3 dias, a economia é de \~1 %, e o ganho é operacional.
   - Empata com o controle autônomo convencional, com \~25 % menos operações por amostra.
2. **Massa de asteroides em simulação.**
   - M = ⟨h⟩²/(G·p\_geom) iguala a OD em lote até \~1 % em Bennu, com 300 a 2 000× menos processamento.
   - Em Psyche, a correção fechada de J2 dá 0,001 % na órbita A. Durante a espiral A→B1, GM < 0,03 % e J2 \~10 % a bordo a partir do 4º dia. Com o ROOP integrado e sem massa prévia, a faixa de densidade 1–8 g/cm³ se fecha em ±0,2 % no primeiro dia, e a nave chega à Orbit B1 no dia 18 em vez do dia 24–25, com o mesmo Δv; a massa do modelo de forma, sem correção, a deixaria 60 km acima do alvo. Nas pernas B1→D e D→C o ganho se repete (7,3 e 7,0 dias em média), somando \~20 dias na missão; abaixo da B1, a descida precisa ser retrógrada para cruzar as ressonâncias 3:2 e 1:1 de rotação.
3. **Massa de asteroides com dados reais.**
   - Na trajetória reconstruída da OSIRIS-REx, o estimador dá GM = 4,897 m³/s² com uma órbita, contra 4,892 ± 0,006 m³/s² (+0,11 %).
   - Com 16 h de dados, o erro fica entre −0,05 % e +0,30 %, contra −0,36 % da 3ª lei de Kepler com 61 h.
   - O mesmo ajuste mediu o deslocamento da órbita por pressão de radiação (78 m, \~6 × 10⁻⁸ m/s²).
4. **Valor.**
   - Na missão Psyche, 2 a 24 dias realocáveis para a órbita de mapeamento elementar (+1 a +15 %).
   - Numa expedição de mineração, 3 a 15 dias a mais de operação, ou US$ 0,6 a 45 milhões de receita bruta sob as premissas de Hein et al.
5. **Próximos passos.**
   - Aplicar o estimador aos dados brutos de navegação da OSIRIS-REx (PDS), para uma medição independente.
   - Estender a correção fechada ao C22 em rotação.
   - Testar ΔF em GEO e em constelações.
   - Fazer a revisão bibliográfica de anterioridade das contribuições 3 a 5.

## Referências, fontes e código

**Fontes consultadas**

*Base teórica*

- Maia, K. M. (2026). [The Orbital Coherence Function: A Keplerian Pre-Filter for Ensemble Orbit Determination](https://zenodo.org/records/20248439), v2, Zenodo, doi:10.5281/zenodo.20248439, seção 2.2.

*Dados de voo*

- NASA/JPL. [Horizons System](https://ssd.jpl.nasa.gov/horizons/), vetores de estado da OSIRIS-REx (−64) relativos a Bennu (@2101955), eclíptica J2000, 3–7/01/2019, cadência de 2 h; consultado em 21/11/2026.
- Scheeres, D. J. et al. (2019). [The dynamic geophysical environment of (101955) Bennu based on OSIRIS-REx measurements](https://www.boulder.swri.edu/~bottke/Reprints/Scheeres_2019_Nature_Astro_3_352_Dynamic_Environ_Bennu.pdf), Nature Astronomy 3, 352–361 (GM = 4,892 ± 0,006 m³/s²).
- Barnouin, O. et al. (2019). [Main shape parameters of Bennu](https://meetingorganizer.copernicus.org/EPSC-DPS2019/EPSC-DPS2019-140-1.pdf), EPSC-DPS 2019 (elipsoide ajustado 252,78 × 246,20 × 228,69 m).
- Leonard, J. M. et al. [OSIRIS-REx Orbit Determination Performance During the Navigation Campaign](https://ntrs.nasa.gov/api/citations/20190029121/downloads/20190029121.pdf), AAS 19-714.
- NASA/OSIRIS-REx. [Preliminary Survey Phase Plan Narrative](https://sbnarchive.psi.edu/pds4/orex/orex.mission_v8.0/document/mission_phase_plans/Preliminary_Survey_Phase_Plan_Narrative.pdf).
- JAXA/ISAS. [Operação de descida para medida de gravidade da Hayabusa2](https://www.isas.jaxa.jp/en/topics/001799.html).

*Psyche*

- [Mass, Density, and Radius of Asteroid (16) Psyche from High-precision Astrometry](https://doaj.org/article/5f5f7cfa65654311bb5c9f4f3b0e5f04), The Astronomical Journal, 2024 (GM = 1,601 ± 0,017 km³/s²).
- [Dynamical Environment and Surface Characteristics of Asteroid (16) Psyche](https://arxiv.org/pdf/2003.10310) (momentos de inércia e C20, C22; período de rotação 4,195948 h).
- [Psyche Mission Description and Design Rationale](https://pmc.ncbi.nlm.nih.gov/articles/PMC12521304/) (durações das órbitas e transferências).
- [Psyche science operations](https://link.springer.com/article/10.1007/s11214-025-01218-x), Space Science Reviews, 2025 (portões de conhecimento de gravidade).
- Zuber, M. T. et al. (2022). [The Psyche Gravity Investigation](https://link.springer.com/article/10.1007/s11214-022-00905-3), Space Science Reviews.
- FU Berlin. [Psyche mission instruments and orbits](https://www.geo.fu-berlin.de/en/geol/fachrichtungen/planet/projects/psyche/_content/text2.html) (altitudes A, B, C, D).

*Controle orbital e economia*

- Lamy, A. et al. [Autonomous Orbit Control experiment on DEMETER](https://issfd.org/ISSFD_2009/OperationsII/Lamy.pdf), ISSFD 2009.
- De Florio, S.; D'Amico, S.; Radice, G. [Flight Results of Precise Autonomous Orbit Keeping Experiment on PRISMA Mission](https://elib.dlr.de/80873), JSR, 2013.
- Hein, A. M.; Matheson, R.; Fries, D. [A Techno-Economic Analysis of Asteroid Mining](https://www.arxiv.org/pdf/1810.03836), Acta Astronautica 168, 2020.
- [Preço da platina, 19/11/2026](https://metalcharts.org/platinum-price-per-ounce) (US$ 43 666/kg).

*Bibliografia de base*

- Vallado, D. A. *Fundamentals of Astrodynamics and Applications*, 4ª ed. Microcosm Press, 2013.
- Goldstein, H.; Poole, C.; Safko, J. *Classical Mechanics*, 3ª ed. Addison-Wesley, 2002.

*Código e dados (Python, sementes fixas)*

- `orex_bennu_vectors.csv`, `orex_real.py` e `orex_kepler3.py`: validação com a trajetória real da OSIRIS-REx (seção 5.6).
- `sim_bodies.py`, `psyche_j2.py`, `psyche_spiral.py`, `psyche_time.py`, `psyche_roop.py` e `psyche_roop2.py`: Bennu e Psyche com campo de grau 2 em rotação, espiral, linha do tempo e ROOP integrado nas três transferências (seções 5.4, 5.5, 5.7 e 5.8).
- `sim_mass2.py` e `od_cost.py`: estimador contra OD em lote, órbita e sobrevoo (seção 5.2).
- `sim_market.py`, `econ.py` e `bench.py`: Objetivo 1 e economia (seção 4).
- `sim_leo.py`, `sim_star.py`, `sim_mass.py` e `mass_j2.py`: testes negativos (seção 6).
