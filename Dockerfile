FROM ubuntu:jammy
RUN apt update && \
    apt install -y build-essential cmake gdb libopenmpi-dev git xterm && \
    useradd mpi-user && \
    mkdir /app && \
    chown -R mpi-user /app
USER mpi-user
WORKDIR /app

# Run Docker with xterm activated:
# docker container run --rm -e DISPLAY=host.docker.internal:0 -v ${PWD}:/app -it mpi:latest