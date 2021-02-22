[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)

For future me: All scripts in here that do not start with a _ and end with a `.py` file extension will be loaded. Just make sure that they include the decorators provided by `cli`:

```
from edmund.entrypoint import cli

@cli.command()
def some_command():
    ...
```
