function a = exist(s)

// exist - test if a variable or a file exist
//
//  a = exist(s);
//
//  Copyright (c) 2008 Gabriel Peyre

a = exists(s);

if a==0
    [fid,err] = mtlb_fopen(s, 'r'); // mtlb_fopen
    if fid>0
        a = 2;
        mclose(s);
    end
end

endfunction