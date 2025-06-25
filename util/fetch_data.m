function fetch_data(remote_dir, opts)
arguments
    remote_dir string
    opts.local_dir string = string(pwd)+"/";
end
    lstring = "pscp -r" + " """+remote_dir+""""+" """+opts.local_dir+"""";
    system(lstring);
end