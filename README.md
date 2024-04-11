# Writing Samples

Welcome! I'm Emily, and this is my collection of writing samples for technical writing. The readme lists my samples and includes the tools I used for them.

# Technical Writing

This section contains my person technical writing projects.

## Heat Equation Project Documentation

For this project, I solved the 2D steady-state heat equation using the finite element and finite difference methods. I documented the process for setting up the environment in WSL to run the code, and I also displayed my code and results. Check out the documention [here](HeatEquation.md).

Tools:
- Vim and Markdown for authoring
- Jupyter Notebook for running the code

# Academic Writing

This section contains my writing from my graduate studies in computational engineering. 

## ANODE Summary

I wrote a summary for a research paper on augmented neural ordinary differential equations. You can check out the [raw latex file](AcademicWriting/anode_main.tex) and the [pdf](AcademicWriting/anode_main.pdf).

Tools:
- LaTeX for authoring
- pdflatex for compiling 

> [!NOTE]
> For this portfolio, I learned how to compile LaTeX files in the terminal instead of using Overleaf. I ran the following code:

```bash
# Install pdflatex
apt-get update
apt install texlive-latex-base

# Compile with pdflatex
pdflatex anode_main.tex
```
