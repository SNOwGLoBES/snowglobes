CC = g++
CFLAGS = -Wall

of: oscillate_flux
all: oscillate_flux pinched

oscillate_flux:
	$(CC) $(CFLAGS) oscillate_flux.C supernova_mixing.C -o oscillate_flux

pinched:
	$(CC) $(CFLAGS) pinched.C supernova_mixing.C -o pinched

clean:
	rm -f oscillate_flux pinched
