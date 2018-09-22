CC = g++
CFLAGS = -Wall

of: oscillate_flux
all: oscillate_flux pinched

oscillate_flux: oscillate_flux.cc supernova_mixing.cc
	$(CC) $(CFLAGS) oscillate_flux.cc supernova_mixing.cc -o oscillate_flux

pinched: pinched.cc supernova_mixing.cc
	$(CC) $(CFLAGS) pinched.cc supernova_mixing.cc -o pinched

clean:
	rm -f oscillate_flux pinched
