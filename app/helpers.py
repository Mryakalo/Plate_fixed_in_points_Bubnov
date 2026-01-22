def ask_parameter(name, params_type, default):
    while True:
        try:
            parameter = input(f'Insert {name}, default: {default}: ')

            if not parameter:
                return default
            parameter = params_type(parameter)
            break
        except ValueError:
            print(f'Input is not {params_type.__name__}')
    return parameter