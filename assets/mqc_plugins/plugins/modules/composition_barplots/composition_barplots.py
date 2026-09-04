#!/usr/bin/env python
""" MultiQC module to parse output from Composition """
from __future__ import print_function
from collections import OrderedDict
import logging
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule
import seaborn as sns
import pandas as pd
import random
import os

# Initialise the logger
log = logging.getLogger("multiqc")

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Composition Barplots",
            anchor="composition_barplots"
        )
        self.css = {
            "assets/style.css": os.path.join(os.path.dirname(__file__), "..", "style.css")
        }
        self.comp_data = dict()
        taxonomy_level = OrderedDict({
            "level-1.csv":"Kingdom",
            "level-2.csv":"Phylum",
            "level-3.csv":"Class",
            "level-4.csv":"Order",
            "level-5.csv":"Family",
            "level-6.csv":"Genus",
            "level-7.csv":"Species"
        })
        data_dict = OrderedDict([(k, None) for k in taxonomy_level.keys()])

        html = """<style> li:hover {cursor: pointer;}
                </style>"""
        
        # Read input files
        for f in self.find_log_files("composition_barplots", filehandles=True):
            if f['fn'] in data_dict:
                df = pd.read_csv(f['f'], index_col=0)
                data_dict[f['fn']] = df
            else:
                log.warning("{} not recogized by composition barplot module".format(f['fn']))

        html_new = ""
        for fn, df in data_dict.items():
            if df is None:
                continue
            # Check if the input is percentage or fraction, if fraction convert to percent
            # The column sum is not perfect 1 or 100 because of floating point
            # Therefore we can only guess by its range
            if df.sum().between(0.9,1.1).all():
                log.debug("{} seems to have fractions instead of percentages. Converting to percentages.".format(fn))
                df = df.mul(100)

            min_values = [10, 5, 0.1, 0]
            min_interface = {"10":"Grouping taxa < 10% into 'Others'",
                             "5": "Grouping taxa < 5% into 'Others'",
                             "0.1": "Grouping taxa < 0.1% into 'Others'",
                             "0": "No grouping"}

            level_id = taxonomy_level[fn]
            level_interface_id = level_id.lower().capitalize()

            button_before = f"""
                        <div class="btn-group">
                            <a type="dropdown"  class="btn btn-secondary dropdown-toggle" data-toggle="dropdown" aria-expanded="false">
                                {level_id} ▼
                            </a>
                            <ul class="dropdown-menu dropdown-menu-right">

                        """
            button = ""
            for minval in min_values:
                # Group taxa below minval
                df_new = df.apply(self.group_taxa, minval=minval).fillna(0)
                # Sort rows by row sum, move 'Others' to the end
                df_new = df_new.loc[df_new.sum(axis=1).sort_values(ascending=False).index, ]
                df_new = pd.concat([df_new.loc[df_new.index != 'Others',], df_new.loc[df_new.index == 'Others',]])

                level_interface_min = min_interface[str(minval)]
                level_interface_id = level_id + str(minval).replace(".", "-")

                color = ['#a1c9f4','#b9f2f0','#ffb482','#cfcfcf','#fbafe4','#8de5a1','#029e73','#fab0e4','#949494','#ca9161','#d55e00','#d0bbff','#debb9b','#56b4e9','#0173b2','#de8f05','#ece133','#cc78bc','#ff9f9b', '#fffea3']

                #Settings for colors
                if len(df_new) > 20:
                    color_21 = sns.color_palette("pastel", len(df_new)-20).as_hex()
                    random.shuffle(color_21)
                    color = color + color_21
                else:
                    color = color[0:len(df_new)]
                df_color = pd.DataFrame(index=df_new.index)
                df_color['color'] = color

                cats = df_color.to_dict('index')
                data = df_new.to_dict('series')
                self.write_data_file(data, 'composition_barplots')
                if len(data) > 0:
                    if fn in self.comp_data:
                        log.debug("Duplicate sample name found! Overwriting: {}".   format(fn))
                    self.add_data_source(f, fn)
                    self.comp_data[fn] = data
                class_name_after = "content-table active"
                click_id = level_interface_id

                config = {
                'id': 'Composition_barplots' + level_interface_id,
                'use_legend': True,
                # Values are already percentages, so no Counts/Percentages switch is offered.
                # cpswitch_c_active must stay unset: with the switch off, a False value makes the
                # flat-image fallback used above 100 samples emit the plot with display:none.
                'cpswitch': False,                       # Show the 'Counts /   Percentages' switch?
                'xlab': 'Samples',                           # X axis label
                'ylab': 'Percentage',
                'ymin': 0,
                'ymax': 100,
                'labelSize': 8,
                'tt_decimals': 3,
                'tt_suffix': '%',
                'tt_percentages': False,
                'height': 1024
                }

                barplot_html = bargraph.plot(self.comp_data[fn], cats, pconfig=config)
                barplot_html = f"<div class='{class_name_after}' id='{level_interface_id}'>" + barplot_html + "</div>"
                html_new += barplot_html

                button += f"""
                          <li><a type="button" thu-nga="#{level_interface_id}" id="{level_interface_id}" onclick=handleClick(event) >{level_interface_min} </a> </li>
                      """

            html = html + button_before + button + "</ul> </div>" 
        script = """
            <script type="text/javascript">
                function funcx()
                    {{
                        document.getElementById('{a}').click();
                  
                    }}
                setTimeout(funcx, 3000)
                function follow(){{
                    contentTables = document.querySelectorAll('.content-table')
                    contentTables.forEach(contentTable => {{
                        if(contentTable.classList.contains('active')){{
                            contentTable.style='display: block';
                        }}else{{
                            contentTable.style='display: none';
                        }}
                    }});
                }}
                function handleClick(e) {{
                    var curr = e.target.getAttribute('thu-nga');
                    var elems = document.querySelectorAll('.content-table');
                    for (var i = 0; i < elems.length; i++) {{
                        if (elems[i].id === curr.substring(1)) {{
                            elems[i].classList.add('active');
                        }} else {{
                            elems[i].classList.remove('active');
                        }}
                    }}
                    follow();
                }};
                follow();
            </script>
        """.format(a=click_id)
        html = html + html_new + script
        if len(self.comp_data) == 0:
            raise UserWarning

        self.add_section(
            name=" ",
            description="Taxa composition plots illustrate the microbial composition at different taxonomy levels from kingdom to species. The interactive figure below shows the microbial composition at species level. Additional composition barplots for other taxonomy levels can be accessed by clicking on the below tabs.",
            plot=html)
        
    def group_taxa(self, col, minval):
        """Helper function to group taxa into Others when below minval for each sample"""
        new_col = col[col>=minval]
        new_col['Others'] = 100 - new_col.sum()
        return new_col.round(3)
