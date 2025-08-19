from setuptools.command.egg_info import egg_info as setuptools_egg_info

class sage_egg_info(setuptools_egg_info):

    def finalize_options(self):
        from sage.env import SAGE_CONF_WHEEL_URL
        if SAGE_CONF_WHEEL_URL:
            self.distribution.install_requires.append(
                f'sage_conf @ {SAGE_CONF_WHEEL_URL}')
        setuptools_egg_info.finalize_options(self)
