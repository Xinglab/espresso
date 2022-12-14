import subprocess
import time


def try_command_once(command):
    completed_process = subprocess.run(command,
                                       check=False,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
    decoded_stdout = completed_process.stdout.decode()
    decoded_stderr = completed_process.stderr.decode()
    if completed_process.returncode != 0:
        return None, 'command:{}\nstdout:\n{}\n\nstderr:\n{}'.format(
            command, decoded_stdout, decoded_stderr)

    return decoded_stdout, None


def try_command(command, retry_interval_seconds):
    errors = list()

    # The final (retry_seconds: None) allows running the command, but
    # without the ability to wait and retry.
    retry_interval_seconds = retry_interval_seconds + [None]
    for retry_seconds in retry_interval_seconds:
        stdout, error = try_command_once(command)
        if not error:
            return stdout, None

        errors.append(error)
        if retry_seconds is None:
            break

        time.sleep(retry_seconds)

    formatted_errors = list()
    for i, error in enumerate(errors):
        formatted_errors.append('attempt {}\n{}'.format(i, error))

    return None, '\n'.join(formatted_errors)
