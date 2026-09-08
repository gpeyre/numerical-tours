# Numerical Tours website

This is the `gh-pages` branch of the Numerical Tours repository. It uses Jekyll and retains the existing public routes and custom domain.

## Build and preview

Use Ruby 3.1 or newer:

```sh
bundle install
bundle exec jekyll build
bundle exec jekyll serve
```

If a system-installed `jekyll` launcher uses a different Ruby version, invoke Jekyll through Bundler using the selected Ruby environment. No Node dependencies are required for the site itself. The old Grunt files are retained as historical theme sources; the current stylesheet is `assets/css/tours.css`.

## Python collection

The home page and `/python/` render the same catalogue from `_data/python_catalog.json`. Search and topic filters run entirely in the browser, with no external search service. Titles, introductions, keywords, notebook pages, downloads, and real figure thumbnails are exported from the content repository by `maintenance/export_website.py` on `master`.

The `python/<tour>/index.html` pages are generated static reading versions. Do not edit their code or prose directly; update and execute the source notebook, then export it again. The first notebook cell supports local and Colab execution. Colab buttons point to `master`, so deploy the tested content before deploying this branch.

The main catalogue presents 58 distinct Python tours. Four historical notebook variants remain in the content branch. Unfinished conversions in `python/todo/` are not advertised as working tours.

## Other languages

`/archive/` links to MATLAB, Julia, and R. Existing language index URLs and MATLAB tour URLs are retained. Legacy pages use a working HTTPS MathJax endpoint and include navigation back to the Python collection. Legacy-language code has not been execution-tested in the Python overhaul.

## Checks

```sh
node maintenance/test-search.cjs
python maintenance/check_site.py _site
```

The search test runs the actual interaction handlers against a DOM fixture. The site check validates generated catalogue pages, links, images, downloads, and Colab targets without requiring a browser. These checks do not replace visual browser testing.
