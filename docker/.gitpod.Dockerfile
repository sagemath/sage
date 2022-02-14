ARG BASE_GITHUB_REPOSITORY=sagemath/sage
ARG BASE_TAG=dev
FROM ghcr.io/${BASE_GITHUB_REPOSITORY}/sage-docker-gitpod-standard-with-targets:${BASE_TAG} as with-targets
RUN sudo rm -rf /var/cache/debconf/* /var/lib/apt/lists/* /tmp/* /var/tmp/*
# Fast doc rebuilds do not work because
# "loading pickled environment... failed; source directory has changed"
# Until this is fixed, we can as well remove the whole documentation, which saves a lot of space.
RUN rm -Rf /home/gitpod/sage-local/share/doc/sage

FROM ghcr.io/${BASE_GITHUB_REPOSITORY}/sage-docker-gitpod-standard-with-system-packages:${BASE_TAG}
COPY --chown=gitpod:gitpod --from=with-targets /home/gitpod/sage/logs  /home/gitpod/sage/logs
COPY --chown=gitpod:gitpod --from=with-targets /home/gitpod/sage-local /home/gitpod/sage-local
