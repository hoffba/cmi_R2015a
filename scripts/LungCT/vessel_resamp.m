function [I,info] = vessel_resamp(I,info,d_new,voxsz_new,interpm)
    
    if ~all(info.voxsz==voxsz_new)

        % Determine new matrix size
        if isempty(d_new)
            d_new = ceil(info.fov./voxsz_new);
        end

        dclass = class(I);
        if ~ismember(dclass,{'single','double'})
            I = single(I);
        end

        % Set up image geometry
        ext = (info.fov - info.voxsz)/2;
        F = griddedInterpolant({linspace(-ext(1),ext(1),info.d(1)),...
                                linspace(-ext(2),ext(2),info.d(2)),...
                                linspace(-ext(3),ext(3),info.d(3))},I,interpm);

        % Interpolate to new voxel locations
        ext = (d_new-1).*voxsz_new/2;
        I = F({linspace(-ext(1),ext(1),d_new(1)),...
               linspace(-ext(2),ext(2),d_new(2)),...
               linspace(-ext(3),ext(3),d_new(3))});
        
        I = cast(I,dclass);

        info.voxsz = voxsz_new;
        info.d = size(I);
    end