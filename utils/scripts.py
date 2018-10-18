

def format_script(file_path, **kwargs):
    formatted_script = file(file_path+".format").read().format(**kwargs)
    file(file_path, "w+").write(formatted_script)