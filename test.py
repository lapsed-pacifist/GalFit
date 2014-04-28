def f():
    while 1:
        try:
            continue
        finally:
            break
print f()