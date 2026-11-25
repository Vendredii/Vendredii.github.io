# Vendredii.github.io

Personal academic homepage, built by Jekyll with no plugins and no
JavaScript, so GitHub Pages builds it automatically on every push.

Pages: Home, News, Projects (one page per project), CV (with publications).

## Where to edit

| To change...                                | Edit                          |
| ------------------------------------------- | ----------------------------- |
| Home page: large photo, name, bio, links    | `_data/profile.yml`           |
| News (the newest one is shown on Home)      | one Markdown file each in `_news/` |
| Projects                                    | one Markdown file each in `_projects/` |
| CV: education, experience, awards           | `_data/profile.yml`           |
| CV: publications                            | `_data/publications.yml`      |
| CV: downloadable PDF                        | put it in `assets/` and set `cv_pdf` in `_data/profile.yml` |
| CV: collaboration network map               | places and people in `_data/collaborators.yml`; links come from the authors in `_data/publications.yml` (script: `assets/js/collab-network.js`) |
| Colours, fonts, page width                  | top of `assets/css/style.css` |
| Header, navigation, footer                  | `_layouts/default.html`       |
| Layout of one page                          | `index.html`, `news.html`, `projects.html`, `cv.html` |

Each data file has comments explaining its fields. A section whose list is
empty or commented out is hidden.

### Adding news

Copy the file in `_news/`, rename it (`YYYY-MM-DD-short-name.md`), and edit
the lines at the top: `title`, `date`, `image` and `summary`. The text below
them is the full article, written in Markdown. The newest item is shown on
the home page; all of them are listed on the News page.

### Adding a project

Copy a file in `_projects/`, rename it, and edit the lines at the top
(`title`, `summary`, `image`, and `order` for its menu position).
The text below those lines is the project page, written in Markdown. The
project appears in the Projects menu and on the Projects page automatically.

### Images

Put images in `assets/images/` and refer to them as
`/assets/images/your-file.jpg`. Other files, such as a CV, can go anywhere
under `assets/`.

## Preview locally

```bash
jekyll serve
```

Then open <http://localhost:4000>. Edits reload automatically, except
`_config.yml`, which needs a server restart.

## Publish

Push to the `Vendredii.github.io` repository. In its Settings > Pages, the
source should be "Deploy from a branch", set to the branch you push to.
