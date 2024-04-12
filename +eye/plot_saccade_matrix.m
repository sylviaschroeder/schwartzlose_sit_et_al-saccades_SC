function plot_saccade_matrix(saccade_matrix)

[nS, nN] = size(saccade_matrix);

imagesc(1:nN, 1:nS, saccade_matrix); 
colormap(GreenWhiteMagenta)
ylabel('Saccade ID')
xlabel('Neuron ID')
title ('In screen saccades')
axis image
formatAxes

end
