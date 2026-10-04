"""Development entry intentionally locked before any observation loading."""
DEVELOPMENT_RELEASED = False

def main():
    if not DEVELOPMENT_RELEASED:
        raise RuntimeError('Development has not been released by the research lead')
    raise RuntimeError('No numerical development driver is installed at PRE-RUN')

if __name__ == '__main__':
    main()
