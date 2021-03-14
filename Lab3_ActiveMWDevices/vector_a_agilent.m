
function vector_a_agilent(senal,nombredatos)

    normaliza=max(max(abs(real(senal))),max(abs(imag(senal))));
    disp(['La señal se ha dividido entre ' num2str(normaliza)]);
    disp('La amplitud del generador es ');
    disp([ num2str(20*log10(normaliza)) ' dB mayor que la potencia real transmitida']);
    senal=senal/normaliza;

    I=real(senal);
    Q=imag(senal);

    sampi=round(I*32768);
    sampq=round(Q*32768);
%     sampi=round(I*32765);
%     sampq=round(Q*32765);

    cuales=find(sampi>32767);
    sampi(cuales)=32767;
%     sampi(cuales)=32750;
    cuales=find(sampq>32767);
    sampq(cuales)=32767;
%     sampq(cuales)=32750;

    % Se escriben en un fichero auxiliar los byte correspondientes a la señal
    fd=fopen(nombredatos,'wb');
    %fd=fopen('tono.wb')
    for contador=1:length(sampi)
      % Saco el MSB de I
      [msbi,lsbi]=parte(sampi(contador));
      % Saco el MSB de Q
      [msbq,lsbq]=parte(sampq(contador));
      fwrite(fd,[msbi lsbi msbq lsbq],'uchar');
    end;
    fclose(fd);

end

function [msb, lsb]=parte(valor)

    if valor>32767
      valor=32767;
    end
    if valor<-32768
      valor=-32768;
    end;
    if valor<0
      valor=valor+65536;
    end;
    msb=bitand(valor,15*4096+15*256)/256;
    lsb=bitand(valor,255);

end
    
