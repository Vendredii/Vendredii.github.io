# Vendredii.github.io

Personal academic homepage. Three pages (Home, Publications, CV) built by Jekyll, with no plugins
and no JavaScript, so GitHub Pages builds it automatically on every push.

## Where to edit

| To change...                                   | Edit                       |
| ---------------------------------------------- | -------------------------- |
| Name, bio, links (home page)                   | `_data/profile.yml`        |
| News (home page shows the newest few)          | `_data/news.yml`           |
| Publications page                              | `_data/publications.yml`   |
| CV page: education, experience, awards, PDF    | `_data/profile.yml`        |
| Colours, fonts, page width                     | top of `assets/css/style.css` |
| Header, navigation, footer                     | `_layouts/default.html`    |
| Layout of one page                             | `index.html`, `publications.html`, `cv.html` |

Each data file has comments explaining its fields. A section whose list is
empty or commented out is hidden.

Images go in `assets/images/`. Other files, such as a CV, can go anywhere
under `assets/` and be linked as `/assets/cv.pdf`.

## Preview locally

```bash
jekyll serve
```

Then open <http://localhost:4000>. Edits to data files, `index.html` and the
stylesheet reload automatically; `_config.yml` needs a server restart.

## Publish

Push to the default branch of the `Vendredii.github.io` repository. In the
repository's Settings > Pages, the source should be "Deploy from a branch".
