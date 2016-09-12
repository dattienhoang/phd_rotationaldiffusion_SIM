PRO tau_log_normal_20130513, mean, FWHM, n_molecules, taus_dist
;Program produces a log normal distribution with a given mean in log value, FWHM for a specific number
;of molecules and returns the taus in the array tau_dist.
;n_molecules=10000
;FWHM = 0.5
;FWHM = 1.2
mean = 2.0
sigma = FWHM/(2.0*(sqrt(2*alog(2))))
;mean = 2
;Gaussian distribution of log(tau)
x=randomN(seed,n_molecules)*sigma + mean

;plot histogram
;binsize = (max(x)-min(x))/100
;binsize = 0.05
;histo_data = histogram(X, binsize = binsize)
;bins = findgen(n_elements(histo_data))*binsize+min(X)
;
;plot, bins, histo_data, YRANGE = [MIN(histo_data)-1, MAX(histo_data)+1],$
; /ystyle, psym=10, xrange = [bins[0], bins[n_elements(histo_data)-1]]$
; , /xstyle
 
;produce log normal distribution of tau values
;if X is a random variable with a normal distribution,
;then Y = exp(X) has a log normal distribution
Y=10.0^(X)
;print, 10^mean(Y)
binsize = 10.0^(0.05)
histo_data = histogram(Y, binsize = binsize)
bins = findgen(n_elements(histo_data))*binsize+min(Y)
plot, bins, histo_data, YRANGE = [MIN(histo_data)-1, MAX(histo_data)+1],$
 /ystyle, psym=10, xrange = [bins[0], bins[n_elements(histo_data)-1]]$
 , /xstyle, /xlog
 
taus_dist = y
get_lun, lun 
filename = strcompress('tau_log_normal_data_FWHM' + string(FWHM))
openw, lun, filename
printf, lun, 'Number of molecules:', n_molecules 
printf, lun, 'st dev of log(tau) distribution:', stdev(alog10(taus_dist))
printf, lun, 'FWHM of log(tau) distribution:',  2*sqrt(2*alog(2))*stdev(alog10(taus_dist))
for i=0, n_molecules - 1 DO printf, lun, taus_dist[i]
close, lun
free_Lun, lun
print, 'tau_log_normal finished'
END