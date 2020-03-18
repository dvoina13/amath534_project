function x = supervisor(params)
    if params.choice == 1 % sine
        
        if ~isfield(params,'amplitude')
            params.amplitude = 1;
        end
        
        if ~isfield(params,'frequency')
            params.frequency = 5;
        end
        
        if ~isfield(params,'phaseshift')
            params.phaseshift = 0;
        end
        
        if ~isfield(params,'verticalshift')
            params.verticalshift = 0;
        end
            
        x = params.amplitude*(sin(2*pi*(1:1:params.duration)*params.frequency*params.sampletime/1000-params.phaseshift)-params.verticalshift);
    elseif params.choice == 2 % sawtooth
        if ~isfield(params,'frequency')
            params.frequency = 5;
        end
        
        %x = sawtooth(2*pi*(1:1:params.duration)*params.frequency*params.sampletime/1000);
        x = asin(sin(2*pi*(1:1:params.duration)*params.frequency*params.sampletime/1000));
    elseif params.choice == 3 % sineproduct
        %x = sin(8?t)sin(12?t)
        x = sin(2*pi*(1:1:params.duration)*4*params.sampletime/1000).*sin(2*pi*(1:1:params.duration)*6*params.sampletime/1000);
    elseif params.choice == 4 % noisysineproduct
        %x = sin(8?t)sin(12?t) + 0.05*xi
        params.choice = 3; %% get sine product and add noise
        x = supervisor(params) + 0.05*randn(1,params.duration);
    elseif params.choice == 5 % oscillator1
        %x = 12 sin(8?t) + 16 sin(12?t) + 14 sin(28?t)
        x = 12*(sin(8*pi*(1:1:params.duration)*params.sampletime/1000))+ 16*(sin(12*pi*(1:1:params.duration)*params.sampletime/1000)) + 14*(sin(28*pi*(1:1:params.duration)*params.sampletime/1000));
    elseif params.choice == 6 % oscillator2
        %x =sin(4?t)sin(6?t)sin(14?t)
        x = (sin(4*pi*(1:1:params.duration)*params.sampletime/1000)).*(sin(6*pi*(1:1:params.duration)*params.sampletime/1000)).*(sin(14*pi*(1:1:params.duration)*params.sampletime/1000)); 
    end
    
    if isfield(params,'normalize') && params.normalize
        x = normalize(x,'range');
    end
    
end