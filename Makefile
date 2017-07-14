all: grid.csv  scheduled_order.csv


grid.csv: bayestar.fits.gz
	python grid.py -o $@ $^

scheduled_order.csv: grid.csv
	python emgwfollowup_scheduling.py -o $@  $^ \
	--method optimal \
	--time '2017/6/12 7:0:0'

clean:
	rm grid.csv scheduled_order.csv
