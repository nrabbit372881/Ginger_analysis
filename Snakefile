configfile: "configs/config.json"

# configure docker mounting options
# ==============================================================================
docker_mount_opt = ""
for volum in config["volumes"]:
    docker_mount_opt += "-v %s:%s:%s " % (volume["real"], volume["virtual"], volume["mode"])
    
    if volume["is_workdir"]:
        docker_mount_opt += "-w %s " % volume["virtual"]