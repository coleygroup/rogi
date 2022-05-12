import rogi


def test_name():
    try:
        assert rogi.__name__ == 'rogi'
        version = rogi.__version__
    except Exception as e:
        raise e