% Copyright David Rowe 2009
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Plot resampled ampltiude modelling information from dump files.

function plresample(samname, f)
  
  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);

  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);

  modelq_name = strcat(samname,"_qmodel.txt");
  if (file_in_path(".",modelq_name))
    modelq = load(modelq_name);
  endif

  resample_name = strcat(samname,"_res.txt");
  if (file_in_path(".",resample_name))
    resample = load(resample_name);
  endif

  k = ' ';
  do 
    figure(1);
    clf;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    plot(s);
    axis([1 length(s) -20000 20000]);

    figure(2);
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    plot((1:L)*Wo*4000/pi, 20*log10(Am),";Am;r");
    axis([1 4000 -10 80]);
    hold on;
    plot((0:255)*4000/256, Sw(f,:),";Sw;");
     
    if (file_in_path(".",modelq_name))
      Amq = modelq(f,3:(L+2));
      plot((1:L)*Wo*4000/pi, 20*log10(Amq),";Amq;g" );
      signal = Am * Am';
      noise = (Am-Amq) * (Am-Amq)'; 
      snr1 = 10*log10(signal/noise);
      Am_err_label = sprintf(";Am error SNR %4.2f dB;m",snr1);
      plot((1:L)*Wo*4000/pi, 20*log10(Amq) - 20*log10(Am), Am_err_label);
      %Am(1:4)
      %Amq(1:4)
    endif

    if (file_in_path(".",resample_name))
	Wo_r = resample(f,1);
        L_r = resample(f,2);
        Am_r = resample(f,3:(L_r+2));
        plot((1:L_r)*Wo_r*4000/pi, 20*log10(Am_r(1:L_r)),"+;Amres;c" );
        %Am_r(1:4)
    endif

    hold off;

    % interactive menu

    printf("\rframe: %d  menu: n-next  b-back  p-png  q-quit", f);
    fflush(stdout);
    k = kbhit();
    if (k == 'n')
      f = f + 1;
    endif
    if (k == 'b')
      f = f - 1;
    endif

    % optional print to PNG

    if (k == 'p')
      figure(1);
      pngname = sprintf("%s_%d_sn.png",samname,f);
      print(pngname, '-dpng', "-S500,500")
      pngname = sprintf("%s_%d_sn_large.png",samname,f);
      print(pngname, '-dpng', "-S800,600")

      figure(2);
      pngname = sprintf("%s_%d_sw.png",samname,f);
      print(pngname, '-dpng', "-S500,500")
      pngname = sprintf("%s_%d_sw_large.png",samname,f);
      print(pngname, '-dpng', "-S1200,800")
    endif

  until (k == 'q')
  printf("\n");

endfunction
