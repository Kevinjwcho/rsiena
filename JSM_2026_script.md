# JSM 2026 — Presentation Script
### *Modeling Longitudinal Three-way Network Data*
**Jinwoo Cho & Nynke Niezink — Carnegie Mellon University**

*Target length: ~14 minutes (~1,950 words). Pace ≈ 140 wpm. Timing markers are cumulative.*

---

## Slide 1 — Title  *(0:00–0:30)*

Good morning, and thank you for coming. I'm Jinwoo Cho, and this is joint work with Nynke Niezink at Carnegie Mellon. Today I'll talk about modeling three-way network data over time — that is, not just who is connected to whom, but how different people *perceive* the network, and how those perceptions evolve. We do this by extending the Stochastic Actor-Oriented Model to Cognitive Social Structures.

---

## Slide 2 — Illusion of the "Objective" Network  *(0:30–2:00)*

Let me begin with an assumption built into most network analysis. We usually analyze one objective network for the sample — a single representative structure that all individuals are assumed to share and agree on. And implicitly, we assume actors are aware of the complete network state — that they know the whole picture.

But that's not how people experience their social world. Each person carries their own subjective map of who is tied to whom, and those maps can be systematically wrong. As a simple example: person A might be convinced that B and C are close friends, even if they barely talk in reality.

This idea has a rich history, and three contributions shaped it: Krackhardt, in 1987, introduced the Cognitive Social Structure — the idea that each person holds a unique view of the network; Casciaro later showed empirically that people in the *same* network perceive it differently depending on their personality and position; and Brands reviewed how these subjective maps form and go on to shape real organizational behavior.

---

## Slide 3 — Why Study Network Perception?  *(2:00–2:50)*

Why does this matter? The Thomas Theorem captures it: "If men define situations as real, they are real in their consequences." It's often the perception, not the underlying truth, that drives behavior — someone may disengage because they *believe* they're isolated, whether or not they actually are.

So we present a statistical model for the dynamics of network perception. Two questions motivate us. What cognitive mechanisms drive how perceptions form and evolve? And — the part I find most interesting — do perceptions of *other people's* ties follow a different mechanism than reports of *our own* ties?

---

## Slide 4 — From Standard Networks to Cognitive Social Structures  *(2:50–4:00)*

Here's the data structure. A standard network is an N-by-N matrix per wave: sender j to receiver k. A Cognitive Social Structure — Krackhardt's idea — adds a third dimension, the **perceiver** i, giving an N-by-N-by-N array observed at several time points.

The entry reads: x-i-j-k equals one if perceiver i believes that sender j has a tie to receiver k. On the right, each slice of the cube is one perceiver's entire view of the network. And notice the red row: that's the case where the perceiver is also the sender, i equals j — the person's report about their *own* outgoing ties, their self-report. So self-reports live inside the same cube; they're just the diagonal slices. That's what lets us compare perception and self-report on equal footing.

---

## Slide 5 — Three-way Stochastic Actor-Oriented Model  *(4:00–5:30)*

How do these perceptions evolve? We extend the SAOM. Change happens between waves through a sequence of unobserved **micro-steps**. In each micro-step, one perceiver i reconsiders a single perceived tie — one j-to-k link — and may flip it. The process is Markov: the next change depends only on the current state.

There are two ingredients. First, the **rate parameters**. We allow two separate rates for each perceiver in each wave: one for how often she revisits her *own* self-reported ties, and a different one for how often she revisits her *perceptions of others* — because there's no reason those two update at the same speed. Once she's activated, a sender j is drawn uniformly at random for her to re-evaluate.

Second, given that sender, she updates the tie j-to-k according to a **choice probability** — this multinomial logit. The probability of toggling a tie is proportional to the exponential of an objective function encoding cognitive tendencies, normalized over possible receivers. So the rate governs *when and who*; the choice probability governs *what*.

---

## Slide 6 — Network Statistics: Reciprocity  *(5:30–6:30)*

What goes into that objective function? A weighted sum of network statistics, each a cognitive mechanism. Let me show reciprocity, because it illustrates the perception-versus-self distinction. Blue arrows are perceived ties, red are self-reported.

On the left, **perceived reciprocity**: inside i's mental map, if she sees j tied to k, does she also tend to see k tied back to j? This captures whether people *assume* relationships are mutual — a balance heuristic for filling in what they can't observe. On the right, **self-reported reciprocity**, the i-equals-j case: if i reports a tie to j, does j tend to report one back? Same structural idea — but one lives in perception, the other in self-report, and we estimate them separately.

---

## Slide 7 — Covariate Effects: Perceived Attributes  *(6:30–7:20)*

We can also bring in actor attributes — the orange nodes carry an attribute like gender. Three effects: a **perceived sender effect** — does the perceiver assume actors with that attribute send more ties; a **perceived receiver effect** — does she assume they receive more; and **perceived homophily** — does she infer ties simply because two people are demographically similar. That last one is a cognitive shortcut we can test directly.

---

## Slide 8 — Estimation via Method of Moments  *(7:20–8:20)*

On estimation: the likelihood is intractable, because we never observe the micro-steps between waves — only the snapshots — and integrating over all possible paths is infeasible. So we use simulation-based Method of Moments.

The logic is simple. We compute target statistics from the real data, T. We simulate the model forward and compute the same statistics on the simulated networks. Then we tune the parameters until the simulated statistics match the observed ones on average — the moment condition, expectation of T-hat equals T. In plain terms: we adjust parameters until the networks our model generates look like the real one. This runs on the same Robbins-Monro machinery RSiena already uses.

---

## Slide 9 — Modeling Flexibility: Homogeneity vs. Heterogeneity  *(8:20–9:30)*

A key modeling choice is *whose* cognitive rules we're estimating. Two extremes.

Under **homogeneity**, we assume some cognitive patterns are shared by everyone. We aggregate statistics across all perceivers and estimate a single, constant parameter — for instance, one perceived-reciprocity effect that applies to the whole sample. That's parsimonious and gives us a population-level statement.

At the other end, **perceiver heterogeneity**: different people use different heuristics. We treat each perceiver's statistics as distinct and allow a separate parameter for each individual — so perceived reciprocity could differ from person to person. In our RSiena implementation, this is exactly the `shareParameters` switch: turn it on to share a rule across perceivers, turn it off to let each perceiver have their own. That flexibility lets us ask whether cognitive biases are universal or individual.

---

## Slide 10 — Empirical Application: Undergraduate Advice Network  *(9:30–10:20)*

Let me make this concrete with real data, from Hunter's 2019 study — a six-week undergraduate summer research program. Twenty students, ten men and ten women, surveyed weekly, giving six waves. The relation is advice-seeking: "Who do you go to for help or advice?"

Because it's a Cognitive Social Structure, each wave is a full twenty-by-twenty-by-twenty array — every student reports their perception of the entire advice network. Our goal is to compare three mechanisms — reciprocity, transitivity, and gender homophily — across perception and self-report.

---

## Slide 11 — Self-Reported Advice Network Across Waves  *(10:20–10:55)*

First, the self-reported network — the ties people claim for themselves — across the six waves. Blue is men, red is women. The main takeaway is stability: the structure stays fairly persistent over the six weeks, so people are reasonably consistent about who they turn to.

---

## Slide 12 — Perceived Network: Two Example Perceivers  *(10:55–11:35)*

Now contrast that with perceptions. Here are two perceivers' views of the same network across the six waves. Two things stand out. They're much sparser — people perceive far fewer ties than are actually reported, because you can't observe everyone. And the two perceivers look quite different from each other — the heterogeneity Casciaro pointed to: same network, genuinely different mental maps.

---

## Slide 13 — Network Statistics over Waves  *(11:35–12:30)*

Summarizing those patterns, two systematic differences emerge. On the left, density: perceived networks are consistently *sparser* than self-reported ones — the green line sits above the blue at every wave. People underestimate how many ties exist. On the right, reciprocity: perceived networks show *higher* reciprocity than self-reported ones. When people are uncertain, they default to assuming mutuality — the balance heuristic. These descriptive patterns motivate the formal model.

---

## Slide 14 — Code: Three-way SAOM in RSiena  *(12:30–13:10)*

Briefly on implementation, since people ask whether this is usable. We've built the three-way model into RSiena. You declare the data as a "threeway" dependent variable — the only new argument — and `shareParameters` controls the homogeneity choice I just described. You add effects like reciprocity, transitivity, and same-gender homophily, targeting either the shared perception slices or the self-reported matrix by name, and estimate through the standard `siena07` call. So it's the familiar RSiena pipeline; only the data object and a couple of arguments are new.

---

## Slide 15 — Empirical Results: Parameter Estimates  *(13:10–14:00)*

And here are the results — perception on top, self-report below. The headline is reciprocity. In perception, the estimate is 1.52 — large and strongly positive. In self-report, it's 0.88 — positive, but much smaller. So people impose far more reciprocity on their *perceptions* of others than actually exists in reported behavior: when they don't know, they assume mutuality. That's a genuine, quantifiable cognitive bias. Gender homophily is positive in both — 0.33 versus 0.29 — so similarity acts as a perceptual shortcut. Transitivity is modest in both.

To conclude: we've built a longitudinal framework for Cognitive Social Structures that detects a real, structured discrepancy between perception and self-report. Going forward we're extending the three-way effects, building pre-estimation diagnostics, and applying the model with Dr. Shihan Li and Dr. Tobias Stark. Thank you — I'd welcome your questions and feedback.

---

*If asked "are cognitive rules shared or individual?" — point back to Slide 9 and the `shareParameters` switch.*
