params.samples_list = ['Sample1', 'Sample2', 'Sample3', 'Sample4']
Channel.from( params.samples_list ).into { samples_list; samples_list2 }

samples_list2.subscribe {
    println "[samples_list2] ${it}"
}

process make_txt {
    tag { "${sample_ID}" }
    executor "local"
    echo true

    input:
    val(sample_ID) from samples_list

    output:
    file "${sample_ID}.txt" into samples_files

    script:

    """
    echo "[make_txt] ${sample_ID}"
    echo "[make_txt] ${sample_ID}" > "${sample_ID}.txt"
    """
}

process make_use_dir {
    echo true

    input:
    file("results_list") from samples_files.collect()

    output:
    file("samples_dir") into samples_dir

    script:
    """

    echo "[make_dir]"
    echo "HERE"${results_list}
    pwd
    mkdir -p samples_dir
    mv ${results_list} samples_dir/
    echo "[use_dir]"
    ls samples_dir
    """
}
