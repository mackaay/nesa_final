---
editor_options: 
  markdown: 
    wrap: 72
---

**Take-home assessment**

The assessment is designed to evaluate how effectively you can do
research and problem-solving using modern tools, while maintaining good
technical practice in coding, version control, and writing.

You may use modern tools, including AI-based tools, search tools, and
coding assistants. However, you are expected to exercise judgement,
verify the correctness of your work, and ensure that your final
submission is reproducible, understandable, and technically sound.

This assessment values both quality and speed. We are interested in how
effectively you can reach a sound, well-justified solution within
realistic time constraints.

## Before the test

Before the test, please prepare a git repository from a previous project
that demonstrates your use of version control in practice. Please ensure
that it is something you are permitted to share and that it does not
contain confidential, proprietary, or sensitive material.

## Instructions

Please create a new git repository for this test before starting Task 1,
and use good version control practice throughout the assessment.

Your submission should include: - a Markdown file containing a link to
the final git repository, with access granted to
`chen.chen@nesa.nsw.edu.au`, and briefly describing the structure of the
submission, how to run any code, and how you consolidated the
repositories; - a single final git repository containing the work
completed during this test; - the pre-existing repository you prepared
before the test, incorporated into the final repository in a way that
preserves the relevant history of both repositories; - the technical
writing sample already prepared before the test, included in the final
repository.

If you use modern tools during the assessment, include a complete log in
the repository describing what tools you used, when you used them, what
prompts or queries you issued where practical, what outputs they
produced in summary. The purpose of this log is not to discourage tool
use, but to help us assess your judgement, workflow, and verification
practices.

**Timely submission** is part of the assessment. Please submit as soon
as you have completed the work.

## Tasks

**1.** Let $X$ be a complex torus of dimension $17$ that contains an
abelian variety of dimension $10$. Recall that the algebraic dimension
of $X$ is the transcendence degree of the field of meromorphic functions
of $X$ over $\mathbb{C}$. What is the smallest possible value of the
algebraic dimension of $X$?

The following paper may be helpful for Task 1:
<https://doi.org/10.48550/arXiv.2106.10308>

**2.** Using `aluminium.csv`, which contains columns `strain` and
`stress`, fit a one-breakpoint Bayesian piecewise linear model to
estimate the yield point of the aluminium sample. Let

$$
y_i \mid \beta_1,\beta_2,\tau,\sigma^2 \sim N(\mu_i,\sigma^2),
$$

where

$$
\mu_i =
\begin{cases}
\beta_1 x_i, & x_i \le \tau,\\
\beta_1 \tau + \beta_2 (x_i-\tau), & x_i > \tau,
\end{cases}
$$

so that the two lines meet continuously at the breakpoint $\tau$ (yield
point). Assume independent priors

$$
\beta_j \mid \sigma^2,\tau \sim N(30000,\,20000^2 \sigma^2), \qquad j=1,2,
$$

a discrete uniform prior on $\tau$ over the candidate breakpoints, and

$$
\sigma^2 \mid \tau \sim \mathrm{Inverse\text{-}Gamma}(4,\,13.12).
$$

Report posterior summaries, including the estimated yield point
(breakpoint) by producing a plot of the observed stress-strain data
together with the fitted posterior mean curve. If time permits, build a
small Shiny app that allows the user to vary the prior settings and
explore how the fitted curve, posterior summaries, and estimated yield
point change.

In your solution, explain how you approached the implementation of the
model and how that approach would change if the dataset were much larger
or if the number of candidate breakpoints increased substantially (in
gigabytes). Discuss the main trade-offs you would consider.

**3.** Consolidate the repository created for this test with the
repository you prepared before the test into a single final repository,
in a way that preserves the relevant history of both repositories.
Briefly explain the workflow you used

**4.** Include your technical writing sample in the submission, together
with a short note explaining its context, your role in writing it, and
the intended audience.

## What we will be looking for

We will assess: - correctness and depth of reasoning; - quality and
clarity of code; - reproducibility and organisation of the submission; -
use of version control in a thoughtful and professional manner; -
quality of written communication; - your ability to use modern tools
effectively while maintaining critical judgement; - your ability to work
efficiently and make strong progress within a limited time.
