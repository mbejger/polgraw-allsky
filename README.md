# Polgraw-allsky documentation source files 

The source files are in the `sources` directory. 

Static site is build with [mkdocs](http://www.mkdocs.org): 

  * `pip install mkdocs` 
  * configuration file is `sources/mkdocs.yml`
  * main commands, in the `sources` directory:
    - `mkdocs serve` serves the website to `http://127.0.0.1:8000` for tests and checks 
    - `mkdocs build` build the sources to `site` (to be copied upstairs to the main repo directory, `cp -r site/* ../`) 


`markdown` will be most probably installed with `mkdocs` in your virtual environment (`pyvenv <venv-name>`). Copy the files from `sources/markdown/extensions` to the path in your installation `... /site-packages/markdown/extensions`. 

Other stuff to install: 

  * `python-markdown-math` for the LaTeX rendering in mathjax

