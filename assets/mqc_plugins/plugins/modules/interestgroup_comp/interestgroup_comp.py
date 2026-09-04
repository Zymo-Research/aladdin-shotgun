#!/usr/bin/env python

"""MultiQC plugin module to plot composition of groups of interest"""

from __future__ import print_function
import logging
from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule
import pandas as pd
import seaborn as sns
import random
from collections import OrderedDict

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Halt execution if we've disabled the plugin
        if config.kwargs.get("disable_plugin", False) is True:
            return None

        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name="Abundances of Groups of Interest", target="", anchor="interestgroup_comp"
        )

        self.data_dict = dict()

        # Find the data files
        for f in self.find_log_files("interestgroup_comp", filehandles=True):
            self.add_data_source(f)
            # Parse the data file
            parsed_data = pd.read_csv(f["f"], sep=",", index_col=0)
            if parsed_data is not None:
                self.data_dict[f["fn"].replace("_groupinterest_comp.csv", "")] = parsed_data
            else:
                log.debug("Could not parse interest group composition data in {}".format(f["fn"]))
                raise UserWarning

        # Nothing found - raise a UserWarning to tell MultiQC
        if len(self.data_dict) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning
        
        log.info("Found {} interestgroup_composition reports".format(len(self.data_dict))) 

        self.data_dict = OrderedDict(sorted(self.data_dict.items()))

        if len(self.data_dict):
            self.write_data_file(self.data_dict, "interestgroup_composition_plots")

        datalist = []

        pconfig = {
            "id": "interestgroup_composition_bargraph",
            "title": "Abundances of Groups of Interest",
            # Values are already percentages, so no Counts/Percentages switch is offered.
            # cpswitch_c_active must stay unset: with the switch off, a False value makes the
            # flat-image fallback used above 100 samples emit the plot with display:none.
            "cpswitch": False,
            "yDecimals": True,
            "ymax": 100,
            "ymin": 0,
            "xlab": "Samples",
            "ylab": "Percentages",
            "tt_decimals": 2,
            "tt_percentages": False,
            "tt_suffix": "%",
            "data_labels" : []
        }

        cats = []

        for category_name, data in self.data_dict.items():
            color = ['#a1c9f4','#b9f2f0','#ffb482','#cfcfcf','#fbafe4','#8de5a1','#029e73','#fab0e4','#949494','#ca9161','#d55e00','#d0bbff','#debb9b','#56b4e9','#0173b2','#de8f05','#ece133','#cc78bc','#ff9f9b', '#fffea3']
            if len(data)>0:
                data_todict = data.to_dict()
                pconfig["data_labels"].append({"name":category_name, "ylab":"Percentages"})
                datalist.append(data_todict)
                #Settings for colors
                if len(data) > 20:
                    color_21 = sns.color_palette("pastel", len(data)-20).as_hex()
                    random.shuffle(color_21)
                    color = color + color_21
                else:
                    color = color[0:len(data)]
                full_datacolor = pd.DataFrame(index=data.index)
                full_datacolor['color'] = color

                cats.append(full_datacolor.to_dict('index'))

        html_content = bargraph.plot(datalist, cats, pconfig) 

        self.add_section(
            name=" ",
            description=("This plot depicts the abundances of groups of taxa of interest defined by a user provided list. "
            "Each tab plots the abundances of specific taxa of interest in each group detected in all samples. "
            "The abundances are in percentages of reads among all microbial reads with assigned taxonomy. "
            "If a group of interest is not listed here, it means no read of that group was detected in any sample."),
            plot = html_content
        )
