
For the most recently compiled PDF version of the manuscript, go to
<https://uwmadison.box.com/s/sbupbphpmbng7wn0eavo1qq883g1ckh4>


## File descriptions:

- `__make.command`: This file compiles the LaTeX document. This is the only way
  users should compile it, since clicking "Compile PDF" in RStudio won't complete all
  the necessary steps and will leave around extraneous files.
  Just click on it (on a mac) to use it.
- `__ms.tex`: The main LaTeX document, mostly containing formatting.
  The only manuscript content in this file are the title, authors, and keywords.
- `##-<Description>.tex` files: These files are for the various sections of 
  the paper. Most of the paper content should be contained in these files.
- `ecology.bst`: This tells LaTeX how we want to format citations.
- `ref.bib`: This contains the reference information.


## Citations in `__ms.tex`:

1. Go to Google Scholar; click the top left, and select Settings; then in the
   Bibliography manager section, check show links to import citations into BibTeX, 
   then Save it. This step only need to be done once.
2. Back to the Google Scholar homepage and search the paper you want to cite. 
   E.g. "Generalized linear mixed models for phylogenetic analyses of community 
   structure", under the item, click Import into BibTeX; then copy the whole new page 
   (Cmd + A then Cmd + C).

```
@article{ives2011generalized,
  title={Generalized linear mixed models for phylogenetic analyses of community structure},
  author={Ives, Anthony R and Helmus, Matthew R},
  journal={Ecological Monographs},
  volume={81},
  number={3},
  pages={511--525},
  year={2011},
  publisher={Wiley Online Library}
}
```

3. Open the `refs.bib` in this folder, and paste the citation information there.

4. To cite the paper in the manuscript, copy its key (e.g. `ives2011generalized` here),
   and the LaTeX cite command (`\citep{ives2011generalized}` here) to where you
   want to cite it. For multiple papers, use `\citep{key1,key2}`.
   To add prefixes or suffixes, use `\citep[e.g.,][, see?]{ives2011generalized}` to
   generate "(e.g., Ives & Helmus 2011, see?)".

