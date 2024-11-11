#!/usr/bin/env python

from setuptools import setup, find_packages


_version = '0.0.1'

setup(
    name='multiqc_custom_plugins',
    version=_version,
    description="Custom MultiQC plugins for shotgun pipeline",
    packages=find_packages(),
    include_package_data=True,
    install_requires=['multiqc==1.14','seaborn==0.12.2','pandas==2.0.1'],
    entry_points={
        'multiqc.templates.v1': [
            'aladdin = plugins.templates.aladdin'
        ],
        'multiqc.modules.v1': [
            'plot_heatmap = plugins.modules.plot_heatmap:MultiqcModule',
            'composition_barplots = plugins.modules.composition_barplots:MultiqcModule',
            'interestgroup_comp = plugins.modules.interestgroup_comp:MultiqcModule'
        ],
        'multiqc.hooks.v1': [
            'before_config = plugins.utils.hooks:before_config',
            'execution_start = plugins.utils.hooks:execution_start'
        ],
    }
)
