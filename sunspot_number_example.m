
%--------------------------------------------------------------------------
%Superposed Epoch Tutorial: Sunspot Number Solar Cycle Example
%Liam Kilcommons, University of Colorado
%5-27-2016
%--------------------------------------------------------------------------
sunspot_data = importdata('SN_d_tot_V2.0.csv');

time = sunspot_data(:,4);
ssn = sunspot_data(:,5);

missing = ssn == -1;
time = time(~missing);
ssn = ssn(~missing);

%Center times for each event
%---------------------------

solar_cycle_max_times = [1829+11/12,1837+3/12,1848+2/12,1860+2/12,1870+8/12,1883+12/12,...
1894+1/12,1906+2/12,1917+8/12,1928+4/12,1937+4/12,1947+5/12,...
1958+3/12,1958+3/12,1968+11/12,1979+12/12,1989+7/12,2000+3/12,2014+4/12];

%Since we don't know if the actual max was nearer the beginning or the end of the month,
%subtract half a month to put the center time in the middle of the specified month
solar_cycle_max_times = solar_cycle_max_times-1/24;

%Time Binning
%------------


n_days_per_bin = 30; %Approx month long bins
nbins = floor(365.25*8/n_days_per_bin);
time_bin_edges = linspace(-4,4,nbins+1);

n_solar_cycles = length(solar_cycle_max_times);

%Preallocate an array in which we will store our final results
%each row will represent a solar cycle, each column, a time bin
ssn_by_event = zeros(n_solar_cycles,nbins);

%Loop through events, and bin sunspot number data into time bins, and assign into
%the results array

for m = 1:n_solar_cycles
  this_solar_cycle_max_time = solar_cycle_max_times(m); 
  
  %Create variables for sunspot number observations for just this solar cycle, with 0 epoch time being
  %the sunspot maximum
  this_solar_cycle_times = time > this_solar_cycle_max_time-4 & time < this_solar_cycle_max_time+4;
  this_epoch_time = time(this_solar_cycle_times)-this_solar_cycle_max_time;
  this_ssn = ssn(this_solar_cycle_times);
  
  for t = 1:nbins-1
    in_bin_t = this_epoch_time >= time_bin_edges(t) & this_epoch_time < time_bin_edges(t+1);
    ssn_by_event(m,t) = mean(this_ssn(in_bin_t));
  end
  
end

%Take the median of each column, to represent the average response for that
%time bin, also take the 25th and 75th percentiles of each column, to
%represent the spread in sunspot numbers for that time bin
time_bin_centers = (time_bin_edges(1:end-1)+time_bin_edges(2:end))/2;
average_ssn = nanmedian(ssn_by_event,1); %The 1 means to take the median for each column, not each row
first_quartile_ssn = quantile(ssn_by_event,.25,1);
third_quartile_ssn = quantile(ssn_by_event,.75,1);

%Plot the results
f = figure();
a = axes(); hold(a,'on');
plot(a,time_bin_centers,average_ssn,'b-');
plot(a,time_bin_centers,first_quartile_ssn,'r--');
plot(a,time_bin_centers,third_quartile_ssn,'r--');
xlabel('Years since maximum sunspot number');
ylabel('International Sunspot Number');
title(sprintf('Superposed Epoch Analysis of %d Solar Cycles (%d-%d)',n_solar_cycles,...
    floor(time(1)),floor(time(end))));
print(f,'-dpng','sunspot_number_sea.png');



