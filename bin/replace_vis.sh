#!/usr/bin/env bash
read -r -a HTML_LIST <<< ${1}
HTML_NUM=$(find . -name "*.html" | wc -l)
HTML_NUM=$(( $HTML_NUM-1 ))
for i in $(seq 0 $HTML_NUM);
do

    EMPTYFILE=0
    PREFIX=$(echo ${HTML_LIST[${i}]} | sed 's/-ancombc-barplot.html//')
    REFGROUP=$(cat refgroup.txt)
    [ ! -s "${PREFIX}-ancombc-barplot.html" ] && EMPTYFILE=1

    cp ${HTML_LIST[${i}]} "${PREFIX}_mqc.html"

    if [ $EMPTYFILE -lt 1 ]; then
        sed -i "s/<div id=\"vis\">/<div id=\"vis_${i}\">/" "${PREFIX}_mqc.html"
        sed -i "s/const el = document.getElementById('vis');/const el = document.getElementById('vis_${i}');/" "${PREFIX}_mqc.html"
        sed -i "s/vegaEmbed(\"\#vis\", spec, embedOpt)/vegaEmbed(\"\#vis_${i}\", spec, embedOpt)/" "${PREFIX}_mqc.html"

        sed -i '1i <!--' "${PREFIX}_mqc.html"
        sed -i "2i id: ${PREFIX}" "${PREFIX}_mqc.html"
        sed -i '3i parent_id: ANCOM_Results' "${PREFIX}_mqc.html"
        sed -i '4i parent_name: "ANCOM Differential Taxa Diversity Results"' "${PREFIX}_mqc.html"
        sed -i "5i parent_description: '<a href=\"https://pubmed.ncbi.nlm.nih.gov/26028277/\">ANCOM-BC</a> is a software designed to find differentially expressed taxa between user specified groups. By default, the first alphabetical group is used as the reference in the below plots. Positive log fold change values (<b>enriched</b> in the legend) indicate that the taxa is expressed more commonly in the subject group as opposed to the reference group. Negative log fold change values (<b>depleted</b> in the legend) indicate that the taxa is expressed more commonly in the reference group as opposed to the subject group.  Black lines in the center of each bar depict the standard error range of log fold changes.'" "${PREFIX}_mqc.html"
        sed -i "6i section_name: ${PREFIX} vs ${REFGROUP} Results" "${PREFIX}_mqc.html"
        sed -i "7i description: ANCOM-BC differential taxa diversity results for subject group '${PREFIX}' compared to the reference group '${REFGROUP}'." "${PREFIX}_mqc.html"
        sed -i '8i -->' "${PREFIX}_mqc.html"
    else
	echo '<!--' > "${PREFIX}_mqc.html"
        echo "id: ${PREFIX}" >> "${PREFIX}_mqc.html"
        echo 'parent_id: ANCOM_Results' >> "${PREFIX}_mqc.html"
        echo 'parent_name: "ANCOM Differential Taxa Diversity Results"' >> "${PREFIX}_mqc.html"
	echo "parent_description: '<a href=\"https://pubmed.ncbi.nlm.nih.gov/26028277/\">ANCOM-BC</a> is a software designed to find differentially expressed taxa between user specified groups. By default, the first alphabetical group is used as the reference in the below plots. Positive log fold change values (<b>enriched</b> in the legend) indicate that the taxa is expressed more commonly in the subject group as opposed to the reference group. Negative log fold change values (<b>depleted</b> in the legend) indicate that the taxa is expressed more commonly in the reference group as opposed to the subject group.  Black lines in the center of each bar depict the standard error range of log fold changes.'" >> "${PREFIX}_mqc.html"
        echo "section_name: ${PREFIX} vs ${REFGROUP} Results" >> "${PREFIX}_mqc.html"
        echo "description: No differential taxa diversity results for group '${PREFIX}' compared to reference group '${REFGROUP}' met the significance cutoff." >> "${PREFIX}_mqc.html"
        echo '-->' >> "${PREFIX}_mqc.html"
    fi

done
