//Find max available taxonomy assignement levels
//Adjust the requested taxonomy levels for exporting accordingly
process FIND_MAX_AVAILABLE_TAX{

    input:
    path(ch_tax_tsv)

    output:
    stdout emit: max_tax
    
    script:
    """
    cut -f2 $ch_tax_tsv | awk -F';' '{print NF}' | sort -nu | tail -n 1
    """
}
