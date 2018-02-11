using WAV

function load_sound(file, n0)

    x, fs = wavread(file)
    n = length(x)

    if contains(file,"bird.wav")
        sel = [1:6000; 12500:15000; 22500:24000; 32500:34000]
        keep = filter(x -> !(x in sel), 1:length(x))
        x=x[keep]
    end

    if (n0 != 0) && (n0 < n)
        x = x[1:n0]
    end

    return x/maximum(x)
end
