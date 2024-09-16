#!/bin/bash

CONDA_ENV_FILE="environment.yml"
CONDA_SNAPSHOT_FILE="environment.yml"


conda_env_init () {
    local conda_prefix=$1
    local conf_dir=$2

    rm -Rf "${conda_prefix}"

    if [ "$3" == "--exp-solver" ]; then
        conda env create --prefix "${conda_prefix}" --file "${conf_dir}/${CONDA_ENV_FILE}"  \
                         --force --experimental-solver=libmamba
    else
        conda env create --prefix "${conda_prefix}" --file "${conf_dir}/${CONDA_ENV_FILE}"  \
                         --force
    fi
}


conda_env_update () {
    local conda_prefix=$1
    local conf_dir=$2

    rm -Rf "${conda_prefix}"
    if [ "$3" == "--exp-solver" ]; then
        conda env update --prefix "${conda_prefix}" --file "${conf_dir}/${CONDA_SNAPSHOT_FILE}" \
                         --prune --experimental-solver=libmamba
    else
        conda env update --prefix "${conda_prefix}" --file "${conf_dir}/${CONDA_SNAPSHOT_FILE}" \
                         --prune
    fi
}


conda_generate_snapshot() {

    local conda_prefix=$1
    local conf_dir=$2

    local channels_section_start="^channels:"
    local channels_section_end="^dependencies:"

    local pip_section_start="- pip:"
    local pip_section_end="^prefix:"


    copy_env_section() {
        local section_start=$1
        local section_end=$2

        local temp_file="$conf_dir/.conda_channel_list.tmp"

        ## Copy a section from the original evironment file to the snapshot file.
        sed -n "/$section_start/,/$section_end/p" "$conf_dir/$CONDA_ENV_FILE" \
            | grep -vE '^\s*#|^\s*$' \
            > "$temp_file"
        # sed -i -E "s/\s+-/  -/g" "$temp_file"
        sed -i -e "/$section_start/,/$section_end/!b" \
            -e "/$section_end/!d;r $temp_file" -e 'd' "$conf_dir/$CONDA_SNAPSHOT_FILE"
        rm "$temp_file"
    }


    rm -f "$conf_dir/$CONDA_SNAPSHOT_FILE"
    conda env export --no-builds --prefix "$conda_prefix" --file "$conf_dir/$CONDA_SNAPSHOT_FILE"
    # conda env export --from-history --prefix "$conda_prefix" --file "$conf_dir/$CONDA_SNAPSHOT_FILE"

    ## Add a heading to the snapshot file
    sed -i '/^name:/ { c \
###############################################################################\
# THIS FILE WAS AUTOMATICALLY GENERATED.\
###############################################################################\
\
name: null
            }' "$conf_dir/$CONDA_SNAPSHOT_FILE"

    ## Delete the prefix line
    # sed -i '/^prefix:.*/d' "$conf_dir/$CONDA_SNAPSHOT_FILE"

    ## Copy the channel list from the original evironment file to the snapshot file.
    ## This is to work around the bug that `conda env export`` does not preserve
    ## channel priority.
    copy_env_section "$channels_section_start" "$channels_section_end"

    ## Copy the list of pip packages from the original evironment file to the snapshot file.
    ## This is to work around the bug that `conda env export`` does not preserve Github links.
    copy_env_section "$pip_section_start" "$pip_section_end"

}
