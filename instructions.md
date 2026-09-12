Group Project Instructions for STAT 638 Applied Bayesian Analysis
Purpose: The project provides experience conducting a complete Bayesian analysis of real data.
Each group of 3–4 students will formulate a scientific question, conduct and evaluate a Bayesian
analysis, introduce a meaningful extension, compare the Bayesian analysis with an appropriate
alternative method, and present the results in the context of scientific questions.
The project consists of a written R Markdown report with HTML output and a live presentation
with Q&A.
I. Written Report
Format
Submit:
1. R Markdown source file (.Rmd);
2. rendered HTML file (.html); and
3. any files needed to reproduce the analysis.
The main report is limited to 2,500 words, excluding references, figure/table captions,
supplementary material, and the AI Use Statement.
Code must be included in the .Rmd file so that the analysis is reproducible, but routine code
should be hidden or collapsible in the HTML report or placed outside the 2500 word limit as an
appendix. Avoid raw R output; present results using well-formatted tables, figures, and text.
The report should contain the following sections.
1. Scientific Question and Data
200–250 words
State a clear scientific question. Briefly describe the data, source, outcome, important predictors,
sample size, and any features relevant to the analysis.
2. Bayesian Model
250–300 words
Specify the Bayesian model mathematically, including the likelihood and important parameters.
Define the parameters and explain why the model is appropriate for the scientific question.
3. Prior Distributions
250–300 words
Specify and justify the priors for the important parameters. Explain whether they are informative,
weakly informative, or otherwise motivated. Simply using a package’s default choice is not
sufficient justification.
4. Posterior Computation and Diagnostics
250–300 words
Explain how the posterior distribution was obtained. If MCMC is used, briefly describe the
method and provide appropriate convergence/computational diagnostics. Explain why the
posterior results can be trusted.
5. Posterior Results and Interpretation
350–400 words
Present the main posterior results using appropriate quantities such as posterior means/medians,
credible intervals, posterior probabilities, or predictions. Interpret the results in terms of the
original scientific question.
6. Sensitivity Analysis and Model Checking
250–300 words
Conduct at least one meaningful sensitivity analysis involving the prior or an important modeling
assumption. When appropriate, use posterior predictive or other model checks. Explain whether
the main conclusions are robust.
7. Alternative Analysis, Comparison, and Bayesian Justification
400–450 words
Analyze the same scientific question using an appropriate frequentist or other established
alternative method.
Compare the Bayesian and alternative analyses in terms of relevant estimates, uncertainty,
assumptions, predictions, interpretation, and scientific conclusions.
8. Innovation or Extension
200–250 words
Include at least one meaningful extension beyond a routine application of a standard Bayesian
model. You are not expected to develop a new statistical methodology.
Possible examples include:
• a problem-specific or informative prior;
• a hierarchical extension;
• an alternative likelihood for overdispersion, heavy tails, excess zeros, etc.;
• nonlinear effects or scientifically meaningful interactions;
• Bayesian model averaging/model uncertainty;
• a problem-specific posterior predictive quantity or Bayesian decision rule;
• incorporation of external information;
• Bayesian treatment of missing data; or
• investigation of subgroup heterogeneity.
Explain what motivated the extension, what you changed, and what was learned from it. Simply
using a sophisticated package, machine-learning method, or AI tool does not constitute
innovation.
Finally, answer explicitly:
Why is the Bayesian approach justified for this problem?
Explain what the Bayesian analysis provides beyond the alternative analysis. Do not
automatically conclude that Bayesian analysis is superior. If it provides little practical advantage
for your problem, say so and explain why.
II. Reproducibility and AI Use
The submitted .Rmd file must reproduce the important analyses, tables, and figures in the HTML
report.
Generative AI may be used for coding, debugging, mathematical/statistical assistance, and
writing. However, the group is responsible for the correctness of everything submitted.
Include an AI Use Statement (maximum 50 words) stating:
1. which AI tools were used;
2. how they were used; and
3. how the group verified AI-generated statistical or computational information.
References and supplementary material do not count toward the 2,500-word limit.
III. Project Presentation
Each group will have:
• 8 minutes for presentation
• 5 minutes for Q&A
Use approximately 8–10 substantive slides, excluding the title and reference slides.
A recommended organization is:
1. Scientific Question and Data
2. Bayesian Model
3. Prior Distributions
4. Posterior Computation and Diagnostics
5. Bayesian Results
6. Bayesian Results/Interpretation
7. Sensitivity Analysis and Model Checking
8. Alternative Analysis and Bayesian Comparison
9. Innovation/Extension
10. Conclusions and Justification for the Bayesian Approach
You may combine slides as appropriate.
Slide Guidelines
Slides should emphasize statistical reasoning, results, and interpretation.
• Do not include R code or screenshots of R output.
• Use readable figures, tables, and equations.
• Avoid lengthy mathematical derivations and excessive text.
• Clearly compare the Bayesian and alternative analyses.
• For the innovation slide, clearly communicate:
motivation → extension → what was learned.
• End by answering the scientific question and explaining whether the Bayesian approach
offered a meaningful advantage.
Every group member must participate meaningfully in the presentation.
IV. Questions and Individual Understanding
During the Q&A, questions will be directed to every member of the group and questions can be
from any part of the project. Therefore, every group member is expected to understand the entire
project.