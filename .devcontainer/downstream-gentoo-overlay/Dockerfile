# name the portage image
FROM gentoo/portage:latest as portage

# based on stage3 image
FROM gentoo/stage3:latest

# copy the entire portage volume in
COPY --from=portage /var/db/repos/gentoo /var/db/repos/gentoo

# settings for using binary packages
RUN echo 'FEATURES="${FEATURES} binpkg-request-signature"' >> /etc/portage/make.conf
RUN getuto

# Using default rust-bin instead of rust leads to issues when trying to install python packages interfacing with rust, the rust compiler is not found.
RUN mkdir -p /etc/portage/package.mask
RUN echo "dev-lang/rust-bin" > /etc/portage/package.mask/rust-bin.masked

# continue with image build ...
# set python to 3.11 (and only 3.11)
RUN echo 'PYTHON_TARGETS="python3_11"' >> /etc/portage/make.conf
RUN echo 'PYTHON_SINGLE_TARGET="python3_11"' >> /etc/portage/make.conf
# accept the license for givaro(?)
RUN echo 'ACCEPT_LICENSE="CeCILL-B"' >> /etc/portage/make.conf
# disable portage features not working in containers
RUN echo 'FEATURES="-ipc-sandbox -mount-sandbox -network-sandbox -pid-sandbox"' >> /etc/portage/make.conf
# install the sage-on-gentoo overlay
RUN emerge --noreplace -vg eselect-repository dev-vcs/git
RUN eselect repository enable sage-on-gentoo
RUN emaint sync -r sage-on-gentoo
# install sage-on-gentoo settings
RUN mkdir -p /etc/portage/package.accept_keywords
RUN ln -s /var/db/repos/sage-on-gentoo/package.keywords/sage-devel /etc/portage/package.accept_keywords/sage-devel
RUN mkdir -p /etc/portage/package.use
RUN ln -s /var/db/repos/sage-on-gentoo/package.use/sage /etc/portage/package.use/sage
# set additional "sane" useflags also aiming to minimize some dependencies
RUN ln -s /var/db/repos/sage-on-gentoo/package.use/minimalist.dependencies /etc/portage/package.use/minimalist
# install the sets
RUN mkdir -p /etc/portage/sets
RUN ln -s /var/db/repos/sage-on-gentoo/sets/sagemath /etc/portage/sets/sagemath
RUN ln -s /var/db/repos/sage-on-gentoo/sets/sagemath-optional /etc/portage/sets/sagemath-optional
# install sage
RUN emerge -vg @sagemath
# example install of sage testing a PR
# RUN EGIT_OVERRIDE_REPO_SAGEMATH_SAGE="https://github.com/mkoeppe/sage.git" EGIT_OVERRIDE_BRANCH_SAGEMATH_SAGE="voidlinux_devcontainer" emerge -vg @sagemath
# install optional packages
# RUN emerge -v @sagemath-optional
