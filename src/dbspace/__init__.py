import pkg_resources

my_version = pkg_resources.get_distribution("my-package-name").version

print(f"Using {my_version}")
