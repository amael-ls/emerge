# Project Emerge

## History
This project ran from 2008 to 2013 and aimed to develop new allometric functions for above-ground tree volume and taper volume.

It was initiated by Christine Deleuze (Office National des Forêts).

Unfortunately, not all the objectives of this project were fulfilled. Today, we try to overcome the difficulties encountered during the project and to finally derive a coherent model for different type of volumes.

## Remark on the R code
The code is, so far, in private repo on Gitlab. Depending on the policy of IGN, I will also release it on Github later.

## Creating the document
To compile, it is first necessary to install the following *lua filters*:
- `fontawesome`, for some symbols
```sh
quarto add quarto-ext/fontawesome
```
- `siunitx-quarto`, to translate latex siunitx commands into html
```sh
quarto add amael-ls/quarto_dev/siunitx-quarto
```

Then, just run:
```sh
quarto render
```
in the shell.

## Data
The data are so far not available, but will be released once the project is over
