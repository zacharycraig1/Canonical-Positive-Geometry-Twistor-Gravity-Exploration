FROM sagemath/sagemath:latest

WORKDIR /workspace

# Sage needs python3 present in the container.
RUN sudo apt-get update && sudo apt-get install -y \
    python3 python3-distutils python3-venv python3-pip \
    git tmux htop \
    && sudo apt-get clean

CMD ["bash"]
