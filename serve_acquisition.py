"""Start the public-data development observation sandbox on localhost."""
import argparse,threading
from acquisition_emulator import ObservationService,http_service
from study_joint_review import fixture


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--patient',choices=('CRC02','CRC11','CRC19'),default='CRC19')
    parser.add_argument('--corruption',choices=('random','adversarial'),default='adversarial');args=parser.parse_args()
    service=ObservationService(*fixture(args.patient,args.corruption))
    with http_service(service) as url:
        print(url,flush=True);print('GET /metadata; POST /field or /review with JSON {"id": integer}. Costs and review are simulated.',flush=True)
        try:threading.Event().wait()
        except KeyboardInterrupt:pass


if __name__=='__main__':main()
