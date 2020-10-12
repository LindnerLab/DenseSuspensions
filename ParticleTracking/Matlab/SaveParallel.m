function SaveParallel(path_mat,path_csv,file)
    save(path_mat,'file');
    
    Pall_CSV = [[file.x]', [file.y]', [file.r]'];
    writematrix(Pall_CSV, path_csv,...
        'Delimiter', 'comma');
end