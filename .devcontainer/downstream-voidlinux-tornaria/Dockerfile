FROM ghcr.io/void-linux/void-glibc-full
RUN xbps-install -Suy bash nodejs gnupg2 openssh git
RUN chsh -s /usr/bin/bash
RUN git clone --depth=200 -b sagemath https://github.com/tornaria/void-packages void-sage && cd void-sage && git fetch --depth=500 https://github.com/void-linux/void-packages && git config user.name "Sage Devcontainer" && git config user.email "sage-devcontainer@example.com" && git rebase FETCH_HEAD && export XBPS_ALLOW_CHROOT_BREAKOUT=1 XBPS_CHROOT_CMD=ethereal && ./xbps-src binary-bootstrap && ./xbps-src pkg sagemath && xbps-install --repository=hostdir/binpkgs/sagemath sagemath
