configfile: "configs/config.json"

# configure docker mounting options
# ==============================================================================
docker_mount_opt = ""
for volume in config["volumes"]:
    docker_mount_opt += "-v %s:%s:%s " % (volume["real"], volume["virtual"], volume["mode"])
    
    if volume["is_workdir"]:
        docker_mount_opt += "-w %s " % volume["virtual"] ##what if is_workdir is not true?

# define query function(設計一個尋找用的函數)
# ==============================================================================
def query(**kwargs):
    for d in config["samples"]:
        for k, v in kwargs.items():
            if d[k] == v:
                return d
                
    raise ValueError("No such sample: %s" % kwargs)
    return None

# download references from NCBI dataset
# ==============================================================================
rule ncbi_dload_dehydrated: # dehydrated 模式比較容易 resume 或重跑，可rehydrate
    output:
        fullname="references/ncbi_dataset.zip"
    log:
        "logs/ncbi-dataset/ncbi_dload_dehydrated.log"
    shell:
        # line 37: dirname 是一個 bash 指令，用來從完整路徑中取出「所在的資料夾」，references/ncbi_dataset.zip的references。
        # line 41: 用我的id跑才不會超出權限
        # line 43: 本機先下載（pull）ccc/ncbi-datasets:20230926 這個 Docker image，然後根據這個 image 啟動出一個新的容器，並加上 docker_mount_opt 的 volume 掛載。
        # line 44: 我要docker執行的command
        # line 45 & 46: 1是將 stdout 寫入；2是將 stderr。  寫入兩者等同 &> {log}
        r"""
        mkdir -p $(dirname {output.fullname}) 
        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name ncbi_dload_dehydrated \
            ccc/ncbi-datasets:20230926 \
                datasets download genome accession GCF_018446385.1 \
                    --include gff3,gtf,genome,seq-report \
                    --filename {output.fullname} \
                    --dehydrated \
                2> {log} \
                1> {log} 
        """

rule ncbi_rehydrate:
    input:
        fin="references/ncbi_dataset.zip"
    output:
        dout=directory("references/ncbi_dataset/data")
    log:
        "logs/ncbi-dataset/ncbi_rehydrate.log"
    shell:
        # line 62: 將 {input.fin} 這個 zip 檔解壓縮到它所在的資料夾中；如果檔案已存在則只在需要時更新（-u）；-d 後面接的是目標資料夾路徑，最後要加上要解壓的 zip 檔案路徑。
        r"""
        unzip -u -d $(dirname {input.fin}) {input.fin}

        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name ncbi_rehydrate \
            ccc/ncbi-datasets:20230926 \
                datasets rehydrate \
                    --directory $(dirname {input.fin}) \
                2> {log} \
                1> {log}
        """
               