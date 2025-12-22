FROM gentoo/portage:latest AS portage
FROM gentoo/stage3:latest

# Copy the entire portage volume in
COPY --from=portage /var/db/repos/gentoo /var/db/repos/gentoo

RUN emerge -qv dev-vcs/git dev-lang/python
