# Project Emerge

## History
This project ran from 2008 to 2013 and aimed to develop new allometric functions for above-ground tree volume and taper volume.

It was initiated by Christine Deleuze (Office National des Forêts).

Unfortunately, not all the objectives of this project were fulfilled. Today, we try to overcome the difficulties encountered during the project and to finally derive a coherent model for different type of volumes.

## Code folder
This folder contains the following files:
- 00_notes.qmd

They must be run in the specified order, except `00_notes.qmd` which is an exploratory file intended to explain some of the encountered challenges.

To compile `00_notes.qmd`, it is first necessary to install the following *lua filters*:
- `fontawesome`
```sh
quarto add quarto-ext/fontawesome
```
Then, to get the html, just run:
```sh
quarto render notes.qmd --to html
```
in the shell.

it is also possible to do this directly from R, with the quarto package:
```r
quarto_render(input = "00_notes.qmd", output_format = "html", cache = TRUE, cache_refresh = FALSE)
```
If there is a bug with some R variables, try to set `cache_refresh = TRUE`.

## Data
The original data are stored on the shared host `smb://del1509N015/2024_FairCarbon/data_origin`. Note that from R, maybe the `smb://` is not necessary. The folder `data_origin` is protected and contains the data as provided by Christine Deleuze on 17 September 2024.
