# Sunbelt 2026 — Presentation Script
### *Stochastic Actor-Oriented Modeling for Cognitive Social Structures*
**Jinwoo Cho & Nynke Niezink — Carnegie Mellon University**

*Target length: ~15 minutes (~2,100 words). Pace ≈ 140 wpm. Timing markers are cumulative.*

---

## Slide 1 — Title  *(0:00–0:30)*

Good morning, everyone, and thank you for being here. My name is Jinwoo Cho, and this is joint work with Nynke Niezink at Carnegie Mellon University. Today I want to talk about how we can model not the network itself, but how people *perceive* the network around them — and crucially, how those perceptions change over time. We do this by extending the Stochastic Actor-Oriented Model to what are called Cognitive Social Structures.

---

## Slide 2 — The Illusion of the "Objective" Network  *(0:30–2:00)*

Let me start with an assumption that's built into most of network analysis. When we collect a network, we usually treat it as a single, objective object — one true adjacency matrix that everyone shares. And implicitly, we assume that the actors inside that network are aware of the complete structure: who is tied to whom.

But that's not how people actually experience their social world. In reality, no one has a bird's-eye view of the network. Each person carries their own subjective map of who is connected to whom — and these maps can be systematically wrong.

This idea has a rich history, and I want to highlight three contributions that shaped it. **Krackhardt, in 1987,** introduced the very concept of the Cognitive Social Structure — he was the first to formalize the idea that each individual holds their own unique slice of the overall network, and that we should measure perception, not just ties. **Casciaro, in 1998,** took this empirically and showed that even people embedded in the *same* network perceive it differently, and that those differences track personality and structural position. And **Brands, in 2013,** pulled this literature together in an influential review, showing how these subjective cognitive maps form and, importantly, how they go on to drive real organizational behavior.

So as a simple example: person A might be completely convinced that B and C are close friends — even if, in reality, B and C barely talk. That gap between perception and reality is exactly what we want to model.

---

## Slide 3 — Why Study Network Perception?  *(2:00–3:00)*

Why does this matter? The Thomas Theorem puts it best: "If men define situations as real, they are real in their consequences." In other words, it's often the *perception* that drives behavior, not the underlying truth. An employee may disengage or leave not because they are structurally isolated, but because they *believe* they are.

So in this talk we present a statistical model for the dynamics of network perception. Two questions motivate us. First, what cognitive mechanisms drive how perceptions form and evolve? And second — and this is the part I find most interesting — do perceptions of *other people's* ties follow a different mechanism than reports of *our own* ties? We'll be able to answer that directly.

---

## Slide 4 — From Standard Networks to Cognitive Social Structures  *(3:00–4:15)*

Here's the data structure. A standard network, at each wave, is an N-by-N matrix: sender j to receiver k. A Cognitive Social Structure — again, Krackhardt's idea — adds a third dimension: the **perceiver**, i. So instead of a matrix, each wave is an N-by-N-by-N array, and we observe it at multiple time points.

The entry reads like this: x-i-j-k equals one if perceiver *i* believes that sender *j* has a tie to receiver *k*.

You can see this on the right. Each slice of the cube is one perceiver's entire view of the network. And notice the highlighted red row: that's the special case where the perceiver is also the sender — when i equals j. That row is the person's report about *their own* outgoing ties — their self-report. So self-reports live inside the same cube as perceptions; they're just the diagonal slices. That structure is what lets us compare the two on equal footing.

---

## Slide 5 — Three-way Stochastic Actor-Oriented Model  *(4:15–6:00)*

Now, how do these perceptions evolve? We extend the Stochastic Actor-Oriented Model.

The core idea is that change happens between waves through a sequence of small, unobserved **micro-steps**. In each micro-step, a single perceiver i reconsiders one perceived tie — one j-to-k link — and decides whether to flip it. The process is Markov: the next change depends only on the current state.

There are two ingredients. First, the **rate parameters**. We allow two separate rates for each perceiver i in each wave m: one rate for how often she revisits her *own* self-reported ties, and a different rate for how often she revisits her *perceptions of others*. This separation matters, because there's no reason to assume people update their self-image and their image of others at the same speed. Once a perceiver is activated, a sender j is drawn uniformly at random for her to re-evaluate.

Second, given that sender j, the perceiver updates the tie j-to-k according to a **choice probability** — this multinomial logit form. The probability of toggling a tie is proportional to the exponential of an objective function, which encodes the cognitive tendencies driving the decision. The denominator just normalizes over the possible receivers. So the rate controls *when and who*, and the choice probability controls *what*.

---

## Slide 6 — Network Statistics: Reciprocity  *(6:00–7:15)*

What goes into that objective function? It's a weighted sum of network statistics, each capturing a cognitive mechanism. Let me show one — reciprocity — because it illustrates the perception-versus-self distinction nicely. Throughout, blue arrows are perceived ties and red arrows are self-reported ties.

On the left is **perceived reciprocity**. Inside perceiver i's mental map, if she sees j sending a tie to k, does she also tend to see k sending one back to j? This captures whether people *assume* relationships are mutual — a balance heuristic in how they fill in the gaps of what they can't directly observe.

On the right is **self-reported reciprocity**, the i-equals-j case. Here, if i reports a tie to j, does that tend to come with j reporting a tie back to i? This is reciprocity in the actual reporting behavior. Same structural idea, but one lives in perception and the other in self-report — and we estimate them separately.

---

## Slide 7 — Covariate Effects: Perceived Attributes  *(7:15–8:15)*

We can also bring in actor attributes — here, the orange nodes carry some attribute like gender or race. Three effects matter.

A **perceived sender effect**: does the perceiver assume that actors with a given attribute send more ties? A **perceived receiver effect**: does she assume actors with that attribute receive more ties? And **perceived homophily**: does she infer ties between people simply because they're demographically similar? That last one is a cognitive shortcut — "people like each other tend to be connected" — and it's something we can test directly in the data.

---

## Slide 8 — Estimation via Method of Moments  *(8:15–9:30)*

A quick word on estimation, because the likelihood is intractable here — we never observe the micro-steps between waves, only the snapshots, and integrating over all possible paths is infeasible. So we use simulation-based Method of Moments.

The logic is simple. We pick a set of target statistics computed from the real data — call them T. We then simulate the model forward and compute the same statistics on the simulated networks. We tune our parameters, theta, until the simulated statistics match the observed ones on average — that's the moment condition, expectation of T-hat equals T. In plain terms: we adjust the parameters until the networks our model generates *look like* the real network. This is solved with a Robbins-Monro stochastic approximation, and it's the same machinery RSiena already uses, which matters for our implementation.

---

## Slide 9 — Empirical Application: Undergraduate Advice Network  *(9:30–10:30)*

Let me make this concrete with real data, from Hunter's 2019 study. The setting is a six-week undergraduate summer research program. There are twenty students — ten men and ten women — and the network was collected weekly, giving us six waves. The relation is advice-seeking: "Who do you go to for help or advice?"

Because it's a Cognitive Social Structure, each wave is a full twenty-by-twenty-by-twenty array — every student reports their perception of the *entire* advice network, not just their own ties. Our goal is to study three mechanisms — reciprocity, transitivity, and gender homophily — and to compare how they operate in perception versus self-report.

---

## Slide 10 — Self-Reported Advice Network Across Waves  *(10:30–11:15)*

First, the self-reported network — the ties people claim for themselves — across the six waves. Blue nodes are men, red are women. The main thing to notice is stability: density hovers between about 0.21 and 0.25 the whole time. So the self-reported advice structure is fairly persistent; people are reasonably consistent about who they go to.

---

## Slide 11 — Perceived Network: Two Example Perceivers  *(11:15–12:00)*

Now contrast that with perceptions. Here are two individual perceivers' views of the same network, across the six waves. Two things jump out. First, they're much sparser — people perceive far fewer ties than are actually reported, because you simply can't observe everyone's relationships. And second, the two perceivers look quite different from each other. This is the heterogeneity Casciaro pointed to: same network, genuinely different mental maps.

---

## Slide 12 — Network Statistics over Waves  *(12:00–13:00)*

If we summarize those patterns, two systematic differences emerge. On the left, density: perceived networks are consistently *sparser* than self-reported ones — that's the green line sitting above the blue at every wave. People underestimate the number of ties out there.

On the right, reciprocity: perceived networks show *higher* reciprocity than self-reported networks. That's telling. It suggests that when people are uncertain, they default to assuming relationships are mutual — exactly the balance heuristic I mentioned. They fill in the unknown with "if j helps k, then k probably helps j too." These descriptive patterns motivate the formal model.

---

## Slide 13 — Code: Three-way SAOM in RSiena  *(13:00–13:45)*

A brief note on implementation, since people always ask whether this is usable. We've built the three-way model directly into the RSiena framework. You declare your data as a "threeway" dependent variable — that's the only new argument — and the `shareParameters` option lets you tie effects across perceivers into a single shared cognitive rule. You then add effects, like reciprocity, transitivity, and same-gender homophily, targeting either the shared perception slices or the self-reported matrix by name. Estimation runs through the standard `siena07` call. So the workflow is the familiar RSiena pipeline — only the data object and a couple of arguments are new.

---

## Slide 14 — Empirical Results: Parameter Estimates  *(13:45–14:30)*

And here are the results, with perception on top and self-report on the bottom. Let me focus on the two numbers that matter most.

Look at reciprocity. In perception, the estimate is 1.52 — large and strongly positive. In self-report, it's 0.88 — still positive, but much smaller. That's our headline finding: people impose far more reciprocity on their *perceptions* of others than actually exists in self-reported behavior. When they don't know, they assume mutuality. It's a genuine cognitive bias, and the model quantifies it.

Gender homophily is positive in both — 0.33 in perception, 0.29 in self-report — so demographic similarity acts as a perceptual shortcut: observers are more likely to infer an advice tie between two people of the same gender. Transitivity is positive but modest in both. The overall story is that perception is not just a noisy copy of reality — it's systematically biased toward balance and similarity.

---

## Slide 15 — Conclusion & Future Work  *(14:30–15:00)*

To wrap up. We've established a longitudinal modeling framework for Cognitive Social Structures, and by extending the SAOM to three-way data, we can directly measure the *discrepancy* between how people perceive networks and how they report their own ties — and we find that discrepancy is real and structured.

Going forward, we're extending the catalog of three-way effects, building a pre-estimation diagnostic workflow, and applying the model in collaborative projects. Thank you very much — I'd love your questions and feedback, and my email is on the slide.

---

*Backup slide on homogeneity vs. heterogeneity is available if asked about whether cognitive rules are shared across perceivers or estimated individually.*
