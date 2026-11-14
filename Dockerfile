# L'atelier g2egrid dans un conteneur.
#
#   cd g2egrid
#   docker build -t g2egrid-atelier .
#   docker run --rm -p 127.0.0.1:8765:8765 g2egrid-atelier
#
# Puis ouvrir l'adresse affichee dans les journaux du conteneur : elle porte
# le jeton d'acces (http://127.0.0.1:8765/?jeton=...). Le serveur ecoute sur
# toutes les interfaces DU CONTENEUR ; c'est la publication du port
# (127.0.0.1:8765:8765) qui le garde sur la machine locale.

FROM python:3.11-slim

COPY --from=ghcr.io/astral-sh/uv:0.8 /uv /uvx /bin/

WORKDIR /app
ENV UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1 \
    PYTHONUNBUFFERED=1

# les dependances d'abord : elles changent moins souvent que le code
COPY pyproject.toml uv.lock README.md ./
RUN uv sync --frozen --no-dev --extra api --no-install-project

COPY . .
RUN uv sync --frozen --no-dev --extra api

EXPOSE 8765
CMD ["uv", "run", "--no-sync", "python", "-m", "g2egrid_api", "--hote", "0.0.0.0", "--sans-navigateur"]
